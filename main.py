# -*- coding: utf-8 -*-

import argparse, sys, importlib, time, os, logging
from datetime import datetime
from dask.distributed import Client, as_completed
import settings
import reader
import analysis
import output
import os
import executor
import numpy as np; import pandas as pd
import atoms

logger = logging.getLogger(__name__)


def process(px, **kwargs):

    cube = kwargs.pop('data', '')
    action = kwargs.pop('action', '')
    param = kwargs.pop('param', '')

    pxldrl = atoms.PixelDrill(cube.isel(dict([(param.col_nm, px[1])])).to_series().astype(float), px)

    return action(pxldrl, settings=param)


def main(param):

    start_time = time.process_time()

    cube = reader.ingest(param)

    param.add_dims(cube)
    # param.add_px_list(cube)

    _log_info(logging.getLogger('paramters'), param)
    _log_info(logging.getLogger('cube'), cube)

    if param.col_nm is None and param.row_nm is None and param.dim_nm is not None:
        # single pixel analysis
        import atoms
        import analysis as aa
        import viz

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
        elif localproc and n_workers:
            client = Client(n_workers=n_workers, threads_per_worker=threads_per_worker)
        else:
            client = Client()

        out = output.OutputCointainer(cube, param, name='test_cube')

        s_param = client.scatter(param, broadcast=True)

        try:
            for rowi in range(len(param.row_val)):
                row = cube.isel(dict([(param.row_nm, rowi)])).compute()

                trashold = 250
                q_trashold = 0.5  # TODO change the percentage to be a parameter

                med = np.greater(row.quantile(q_trashold, dim=param.dim_nm), trashold)
                masked = np.ma.MaskedArray(row, np.resize(med, (row.sizes[param.dim_nm], len(med))))
                px_list = list(map(lambda x: [rowi, x], np.argwhere(med.values is False)))

                if px_list:
                    s_row = client.scatter(row, broadcast=True)
                else:
                    continue

                dim_val = pd.to_datetime(
                    param.dim_val).year.unique()  # <-- pd.to_datetime(pd.to_datetime(param.dim_val).year.unique(), format='%Y')
                col_val = range(0, len(param.col_val))
                t_sl = pd.DataFrame(index=dim_val, columns=col_val)
                t_spi = pd.DataFrame(index=dim_val, columns=col_val)
                t_si = pd.DataFrame(index=dim_val, columns=col_val)
                t_cf = pd.DataFrame(index=dim_val, columns=col_val)

                futures = client.map(process, px_list, **{'data': s_row, 'param': s_param, 'action': action})

                for future, pxldrl in as_completed(futures, with_results=True):

                    col = pxldrl.position[1]
                    t_sl.iloc[:, col] = pxldrl.sl[:]
                    t_spi.iloc[:, col] = pxldrl.spi[:]
                    t_si.iloc[:, col] = pxldrl.si[:]
                    t_cf.iloc[:, col] = pxldrl.cf[:]

                    if pxldrl.error:
                        print(pxldrl.position)
                    client.cancel(future)
                    del future, pxldrl, col

                out.sl[rowi] = np.expand_dims(t_sl.transpose().values, axis=0)
                out.spi[rowi] = np.expand_dims(t_spi.transpose().values, axis=0)
                out.si[rowi] = np.expand_dims(t_si.transpose().values, axis=0)
                out.cf[rowi] = np.expand_dims(t_cf.transpose().values, axis=0)

                client.cancel(s_row)
                client.cancel(futures)

                del futures, t_sl, t_cf, t_si, t_spi, row, s_row, px_list
                gc.collect()
            return out

        except Exception as ex:

            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message)

            logger.debug(
                f'Critical error in the main loop, latest position row {rowi}, col {col}, error type {message}')

        result_cube = executor.analyse(cube, client, param, analysis.phenolo, out)

        result_cube.close()
    else:
        raise ValueError

    end = time.process_time() - start_time

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
                                      level=args.log*10,
                                      filemode='w')
        logger = logging.getLogger(__name__)
        logger.info('*** Phenolo 2.0 ***')
        logger.info('Process started @ {}'.format(datetime.now()))

        param = settings.ProjectParameters(path=args.conf, type='ini')  # TODO  type 'ini' must be flexible

        main(param)

    except KeyboardInterrupt:
        print('Process killed')
        raise sys.exit(1)
    except Exception:
        print('Exception occurred')
        raise sys.exit(1)


