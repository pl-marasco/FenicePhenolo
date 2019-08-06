# -*- coding: utf-8 -*-

import argparse
import importlib
import logging
import os
import sys
import time
from datetime import datetime

from dask.distributed import Client

from phenolo import atoms, settings, reader, viz, output, analysis as aa, executor

logger = logging.getLogger(__name__)


def main(param):
    start_time = time.process_time()

    cube = reader.ingest(param)

    _log_info(logging.getLogger('paramters'), param)
    _log_info(logging.getLogger('cube'), cube)

    if param.col_nm is None and param.row_nm is None and param.dim_nm is not None:
        # single pixel analysis

        param.add_px_list(cube)

        if param.col_nm is not None and param.row_nm is not None:
            ts = cube.isel(dict([(param.col_nm, 0), (param.row_nm, 0)])).to_series().astype(float)
        else:
            ts = cube.to_series()

        pxldrl = atoms.PixelDrill(ts, (0, 0))

        if len(param.pixel_list) == 0:
            print('No acceptable value in the position requested')
            sys.exit(1)
        sngpx_pheno = aa.phenolo(pxldrl, settings=param)

        viz.plot(sngpx_pheno)

    elif len(cube.coords.get(param.col_nm)) is not 1 and len(cube.coords.get(param.row_nm)) is not 1:

        localproc = param.processes
        n_workers = param.n_workers
        threads_per_worker = param.threads_per_worker

        if ~localproc and n_workers and threads_per_worker:
            client = Client(processes=localproc, n_workers=n_workers, threads_per_worker=threads_per_worker)
        elif ~localproc and n_workers:
            client = Client(processes=localproc, n_workers=n_workers)
        elif localproc and n_workers:
            client = Client(n_workers=n_workers, threads_per_worker=threads_per_worker)
        else:
            client = Client()

        out = output.OutputCointainer(cube, param, name=param.outName)

        result_cube = executor.analyse(cube, client, param, aa.phenolo, out)

        result_cube.close()

        client.close()
    else:
        raise ValueError

    end = time.process_time() - start_time

    print(f'\rProcess ended @{datetime.now()} after a total time of '
          f'{end} processing {len(param.col_val)*len(param.row_val)}')

    logger.info('Process ended @{} after a total time of {}'.format(datetime.now(), end))


def _log_info(logger, param):
    logger.debug('-------------------- start values --------------------')
    for key, value in param.__dict__.items():
        logger.debug('{} = {}'.format(key, value))
    logger.debug('--------------------  end values  --------------------')


if __name__ == '__main__':

    assert sys.version_info[:2] >= (3, 6), "You need at minimum python 3.6 to execute this script"

    modules = ['numpy', 'pandas', 'xarray', 'rasterio', 'netCDF4', 'scipy', 'pyhdf', 'seasonal', 'dask', 'distributed',
               'matplotlib', ]

    for module in modules:
        assert importlib.util.find_spec(module), "You need {0} module".format(module)

    print('~~~ Phenolo 2.0 ~~~\r\n')

    # Options
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-c', '--conf', help='Configuration file position', required=True)
        parser.add_argument('-p', '--plot', action='store_true', help='Activate the plotting mode', required=False)
        parser.add_argument('-l', '--log', nargs='?', const=2, type=int, help='Activate the log', required=False)
        args = parser.parse_args()

        if (not args.conf) or (not os.path.isfile(os.path.normpath(args.conf))):
            parser.error("the -c argument for the config file is missing or the given path doesn't exist!!!")

        if args.log:
            logpth = os.path.join(os.path.dirname(args.conf), 'log.txt')
            log = logging.basicConfig(filename=logpth,
                                      level=args.log * 10,
                                      filemode='w')
        logger = logging.getLogger(__name__)
        logger.info('*** Phenolo 2.0 ***')
        logger.info('Process started @ {}'.format(datetime.now()))

        param = settings.ProjectParameters(path=args.conf, type='ini')  # TODO  type 'ini' must be flexible

        main(param)

    except KeyboardInterrupt:
        print('Process killed')
        raise sys.exit(1)
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)

        logger.debug(f'Critical error in the main loop, error: {type(ex).__name__, ex.args}')
        print(message)

        raise sys.exit(1)
