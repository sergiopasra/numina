#
# Copyright 2008-2018 Universidad Complutense de Madrid
#
# This file is part of Numina
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""User command line interface of Numina."""

from __future__ import print_function

import datetime
import logging
import os
import errno
import shutil
import pickle

import six
import yaml

import numina.util.objimport as objimp
from numina.types.frame import DataFrameType
from numina.util.context import working_directory


_logger = logging.getLogger(__name__)


class DataManager(object):
    def __init__(self, backend):
        self.backend = backend

    def store_task(self, task, result_dir):

        where = DiskStorageDefault(result_dir)
        where.store(task)


class ProcessingTask(object):
    def __init__(self, obsres=None, insconf=None):

        self.observation = {}
        self.runinfo = {}
        self.result = None
        self.id = 1

        self.time_create = datetime.datetime.utcnow()
        self.time_start = 0
        self.time_end = 0
        self.request = "reduce"
        self.request_params = {}
        self.state = 0

        if insconf:
            self.set_runinfo(insconf)

        self.observation['mode'] = 'unknown'
        self.observation['observing_result'] = 'unknown'
        self.observation['instrument'] = 'unknown'

        if obsres:
            self.set_obsres(obsres)

        if insconf and obsres:
            if insconf['instrument_configuration']:
                self.observation['instrument_configuration'] = insconf['instrument_configuration']

    def set_obsres(self, obsres):
        self.observation['mode'] = obsres.mode
        self.observation['observing_result'] = obsres.id
        self.observation['instrument'] = obsres.instrument

    def set_runinfo(self, insconf):
        self.runinfo['pipeline'] = insconf['pipeline']
        self.runinfo['recipe'] = insconf['recipeclass'].__name__
        self.runinfo['recipe_full_name'] = objimp.fully_qualified_name(insconf['recipeclass'])
        self.runinfo['runner'] = insconf.get('runner', 'numina')
        self.runinfo['runner_version'] = insconf.get('runner_version', "0")
        self.runinfo['data_dir'] = insconf['workenv'].datadir
        self.runinfo['work_dir'] = insconf['workenv'].workdir
        self.runinfo['results_dir'] = insconf['workenv'].resultsdir
        self.runinfo['base_dir'] = insconf['workenv'].basedir
        self.runinfo['recipe_version'] = insconf['recipe_version']

    def store(self, where):

        # save to disk the RecipeResult part and return the file to save it
        saveres = self.result.store_to(where)

        with open(where.result, 'w+') as fd:
            yaml.dump(saveres, fd)

        self.result = where.result

        with open(where.task, 'w+') as fd:
            yaml.dump(self.__dict__, fd)
        return where.task


class WorkEnvironment(object):
    def __init__(self, obsid, basedir, workdir=None,
                 resultsdir=None, datadir=None):

        self.basedir = basedir

        if workdir is None:
            workdir = os.path.join(basedir, 'obsid{}_work'.format(obsid))

        self.workdir = os.path.abspath(workdir)

        if resultsdir is None:
            resultsdir = os.path.join(basedir, 'obsid{}_results'.format(obsid))

        self.resultsdir = os.path.abspath(resultsdir)

        if datadir is None:
            datadir = os.path.join(basedir, 'data')

        self.datadir = os.path.abspath(datadir)
        if six.PY2:
            index_base = "index-2.pkl"
        else:
            index_base = "index.pkl"
        self.index_file = os.path.join(self.workdir, index_base)
        self.hashes = {}

    def sane_work(self):
        # make_sure_path_doesnot_exist(self.workdir)
        _logger.debug('check workdir for working: %r', self.workdir)
        make_sure_path_exists(self.workdir)
        make_sure_file_exists(self.index_file)
        # Load dictionary of hashes

        with open(self.index_file, 'rb') as fd:
            try:
                self.hashes = pickle.load(fd)
            except EOFError:
                self.hashes = {}
        # make_sure_path_doesnot_exist(self.resultsdir)
        make_sure_file_exists(self.index_file)

        # make_sure_path_doesnot_exist(self.resultsdir)
        _logger.debug('check resultsdir to store results %r', self.resultsdir)
        make_sure_path_exists(self.resultsdir)

    def copyfiles(self, obsres, reqs):

        _logger.info('copying files from %r to %r', self.datadir, self.workdir)

        if obsres:
            self.copyfiles_stage1(obsres)

        self.copyfiles_stage2(reqs)

    def copyfiles_stage1(self, obsres):
        import astropy.io.fits as fits
        _logger.debug('copying files from observation result')
        tails = []
        sources = []
        for f in obsres.images:
            if not os.path.isabs(f.filename):
                complete = os.path.abspath(os.path.join(self.datadir, f.filename))
            else:
                complete = f.filename
            head, tail = os.path.split(complete)
            #initial.append(complete)
            tails.append(tail)
#            heads.append(head)
            sources.append(complete)

        dupes = self.check_duplicates(tails)

        for src, obj in zip(sources, obsres.images):
            head, tail = os.path.split(src)
            if tail in dupes:
                # extract UUID
                hdr = fits.getheader(src)
                img_uuid = hdr['UUID']
                root, ext = os.path.splitext(tail)
                key = "{}_{}{}".format(root, img_uuid, ext)

            else:
                key = tail
            dest = os.path.join(self.workdir, key)
            # Update filename in DataFrame
            obj.filename = dest
            self.copy_if_needed(key, src, dest)

        if obsres.results:
            _logger.warning("not copying files in 'results")
        return obsres

    def check_duplicates(self, tails):
        seen = set()
        dupes = set()
        for tail in tails:
            if tail in seen:
                dupes.add(tail)
            else:
                seen.add(tail)
        return dupes

    def copyfiles_stage2(self, reqs):
        _logger.debug('copying files from requirements')
        for _, req in reqs.stored().items():
            if isinstance(req.type, DataFrameType):
                value = getattr(reqs, req.dest)
                if value is None:
                    continue

                complete = os.path.abspath(
                    os.path.join(self.datadir, value.filename)
                )

                self.copy_if_needed(value.filename, complete, self.workdir)

    def copy_if_needed(self, key, src, dest):

        md5hash = compute_md5sum_file(src)
        _logger.debug('compute hash, %s %s %s', key, md5hash, src)

        # Check hash
        hash_in_file = self.hashes.get(key)
        if hash_in_file is None:
            trigger_save = True
            make_copy = True
        elif hash_in_file == md5hash:
            trigger_save = False
            if os.path.isfile(dest):
                make_copy = False
            else:
                make_copy = True
        else:
            trigger_save = True
            make_copy = True
            
        self.hashes[key] = md5hash
            
        if make_copy:
            _logger.debug('copying %r to %r', key, self.workdir)
            shutil.copy(src, dest)
        else:
            _logger.debug('copying %r not needed', key)

        if trigger_save:
            _logger.debug('save hashes')
            with open(self.index_file, 'wb') as fd:
                pickle.dump(self.hashes, fd)

    def adapt_obsres(self, obsres):
        """Adapt obsres after file copy"""

        _logger.debug('adapt observation result for work dir')
        for f in obsres.images:
            # Remove path components
            f.filename = os.path.basename(f.filename)
        return obsres


def compute_md5sum_file(filename):
    import hashlib
    md5 = hashlib.md5()
    with open(filename, 'rb') as f:
        for chunk in iter(lambda: f.read(128 * md5.block_size), b''):
            md5.update(chunk)
    return md5.hexdigest()


def make_sure_path_doesnot_exist(path):
    try:
        shutil.rmtree(path)
    except (OSError, IOError) as exception:
        if exception.errno != errno.ENOENT:
            raise


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except (OSError, IOError) as exception:
        if exception.errno != errno.EEXIST:
            raise


def make_sure_file_exists(path):
    try:
        with open(path, 'a') as fd:
            pass
    except (OSError, IOError) as exception:
        if exception.errno != errno.EEXIST:
            raise


class DiskStorageDefault(object):
    def __init__(self, resultsdir):
        super(DiskStorageDefault, self).__init__()
        self.result = 'result.yaml'
        self.task = 'task.yaml'
        self.resultsdir = resultsdir
        self.idx = 1

    def get_next_basename(self, ext):
        fname = 'product_%03d%s' % (self.idx, ext)
        self.idx = self.idx + 1
        return fname

    def store(self, completed_task):
        """Store the values of the completed task."""

        with working_directory(self.resultsdir):
            _logger.info('storing result')
            return completed_task.store(self)
