import logging


def run():
    log = logging.getLogger('test')

    log.debug('Debug inside test')

    x()

    return


def x():

    logging.debug('z')

    return True
