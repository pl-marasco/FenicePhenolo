# -*- coding: utf-8 -*-
import sys

import gc
import pandas as pd
from dask.distributed import Client, as_completed
import atoms
import logging
import numpy as np
import copy

logger = logging.getLogger(__name__)


class Processor(object):
    def __init__(self):
        pass




def analyse(cube, client, param, action, out):

