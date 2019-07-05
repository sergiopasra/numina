#
# Copyright 2016-2019 Universidad Complutense de Madrid
#
# This file is part of Numina
#
# Numina is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Numina is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Numina.  If not, see <http://www.gnu.org/licenses/>.
#

"""User command line interface of Numina."""


import os
import logging

import six
import numina.drps
from numina.store import load
from numina.dal.absdal import AbsDrpDAL
from numina.exceptions import NoResultFound
from numina.dal.stored import StoredProduct, StoredParameter
from numina.dal.utils import tags_are_valid
from numina.core import DataFrameType
from numina.instrument.assembly import assembly_instrument
from numina.util.fqn import fully_qualified_name

from .db.model import ObservingBlock, DataProduct, RecipeParameters, ObservingBlockAlias
from .db.model import DataProcessingTask, ReductionResult, Frame, ReductionResultValue


_logger = logging.getLogger(__name__)


def search_oblock_from_id(session, obsref):

    # Search possible alias
    alias_res = session.query(ObservingBlockAlias).filter_by(alias=obsref).first()
    if alias_res:
        obsid = alias_res.uuid
    else:
        obsid = obsref

    res = session.query(ObservingBlock).filter(ObservingBlock.id == obsid).one()
    if res:
        return res
    else:
        raise NoResultFound("oblock with id %d not found" % obsid)


DB_FRAME_KEYS = [
    'instrument',
    'object',
    'observation_date',
    'uuid',
    'type',
    'mode',
    'exptime',
    'darktime',
    'insconf',
    'blckuuid',
    'quality_control',
    'vph',
    'insmode'
]


def metadata_fits(obj, drps):

    # First. get instrument
    objl = DataFrameType().convert(obj)

    with objl.open() as hdulist:
        # get instrument
        instrument_id = hdulist[0].header['INSTRUME']

    this_drp = drps.query_by_name(instrument_id)

    datamodel = this_drp.datamodel
    result = DataFrameType(datamodel=datamodel).extract_db_info(obj, DB_FRAME_KEYS)
    return result

def gg():
    class ProcessingTask(object):
        def __init__(self):
            self.result = None
            self.id = 1

            self.time_create = datetime.datetime.utcnow()
            self.time_start = 0
            self.time_end = 0
            self.request = "reduce"
            self.request_params = {}
            self.request_runinfo = self._init_runinfo()
            self.state = 0
            self.obsid = None


        @classmethod
        def _init_runinfo(cls):
            request_runinfo = dict(runner='unknown', runner_version='unknown')
            return request_runinfo
    return 0

def generate_reduction_tasks(session, obid, request_params):
    """Generate reduction tasks."""

    obsres = search_oblock_from_id(session, obid)

    request = {"id": obid}
    request.update(request_params)

    # Generate Main Reduction Task
    print('generate main task')
    dbtask = DataProcessingTask()
    dbtask.host = 'localhost'
    dbtask.label = 'root'
    dbtask.awaited = False
    dbtask.waiting = True
    dbtask.method = 'reduction'
    dbtask.request = request
    dbtask.ob = obsres
    print('generate done')
    session.add(dbtask)
    # Generate reductionOB
    #
    # print('generate recursive')
    # recursive_tasks(dbtask, obsres, request_params)

    session.commit()
    return dbtask

class SqliteDAL(AbsDrpDAL):
    def __init__(self, dialect, session, basedir, datadir, components=None):
        import numina.instrument.assembly as asbl

        drps = numina.drps.get_system_drps()
        com_store = asbl.load_panoply_store(drps)

        super(SqliteDAL, self).__init__(drps)

        self._RESERVED_MODE_NAMES = ['nulo', 'container', 'root', 'raiz']

        self.dialect = dialect
        self.session = session
        self.basedir = basedir
        self.datadir = datadir
        self.components = components if components else com_store
        self.extra_data = {}

    def dump_data(self):
        pass

    def new_task(self, request, request_params):
        # obsres = search_oblock_from_id(session, obid)

        # request = {"id": obid}
        # request.update(request_params)

        # Generate Main Reduction Task
        print('generate main task')
        dbtask = DataProcessingTask()
        dbtask.host = 'localhost'
        dbtask.label = 'root'
        dbtask.awaited = False
        dbtask.waiting = True
        dbtask.request = request
        dbtask.request_params = request_params
        dbtask.request_runinfo = dbtask._init_runinfo()
        # dbtask.ob = obsres
        print('generate done')
        self.session.add(dbtask)
        # Generate reductionOB
        #
        # print('generate recursive')
        # recursive_tasks(dbtask, obsres, request_params)

        self.session.commit()
        return dbtask

    def update_task(self, task):
        _logger.debug('update task=%d in backend', task.id)
        # FIXME: ugly workaround
        # to make the session notice the change, we have to
        # reassign, mutable dictionary
        # check mutabledict
        # newdict_r = task.request_runinfo.copy()
        # newdict_p = task.request_params.copy()
        # newdict = {}
        # task.request_runinfo = newdict
        # task.request_params = newdict
        # self.session.commit()
        # task.request_runinfo = newdict_r
        # task.request_params = newdict_p
        self.session.commit()

    def update_result(self, task, serialized, filename):
        _logger.debug('update result in backend')

        session = self.session
        result = task.result

        DB_PRODUCT_KEYS = [
            'instrument',
            'observation_date',
            'uuid',
            'quality_control'
        ]

        if result is None:
            return

        res_dir = task.request_runinfo['results_dir']

        result_db = ReductionResult()
        result_db.task_id = task.id
        result_db.uuid = str(task.result.uuid)
        result_db.qc = task.result.qc.name
        result_db.mode = task.request_runinfo['mode']
        result_db.pipeline = task.request_runinfo['pipeline']
        result_db.instrument_id = task.request_runinfo['instrument']
        result_db.time_create = task.time_end
        result_db.time_obs = None
        result_db.recipe_class = task.request_runinfo['recipe_class']
        result_db.recipe_fqn =  task.request_runinfo['recipe_fqn']
        result_db.ob_id = task.request_params['oblock_id']
        result_db.result_dir =  res_dir
        result_db.result_file = filename

        session.add(result_db)

        for key, prod in result.stored().items():

            val = ReductionResultValue()

            # fullpath = os.path.join(task.request_runinfo['results_dir'], saveres['values'][prod.dest])
            # relpath = os.path.relpath(fullpath, task.request_runinfo['base_dir'])
            # val.name = prod.dest
            # val.datatype = prod.type.name()
            # val.contents = relpath
            result_db.values.append(val)

            if prod.type.isproduct():

                internal_value = getattr(result, key)
                ometa = prod.type.extract_db_info(internal_value, DB_PRODUCT_KEYS)

                product = DataProduct()
                #product.result_id = result_db.id
                product.qc = task.result.qc.name
                product.instrument_id = task.request_runinfo['instrument']
                product.time_create = task.time_end
                product.time_obs = ometa['observation_date']
                product.uuid = ometa['uuid']
                product.oblock_id = task.request_params['oblock_id']
                product.type = prod.type.name()
                product.type_fqn = fully_qualified_name(prod.type)
                product.content =  os.path.join(res_dir, serialized['values'][key])

                master_tags = ometa['tags']
                for k, v in master_tags.items():
                    product[k] = v

                session.add(product)

        session.commit()

    def add_obs(self, obtable):
        from numina.core.oresult import ObservationResult
        import uuid
        import datetime

        loaded_data = obtable
        obs_blocks = {}
        # complete the blocks...
        for el in loaded_data:
            obs_blocks[el['id']] = el

        # FIXME: id could be UUID
        obs_blocks1 = {}
        for ob_id, block in obs_blocks.items():
            ob = ObservationResult(
                instrument=block['instrument'],
                mode=block['mode']
            )
            ob.id = ob_id
            ob.uuid = str(uuid.uuid4())
            ob.configuration = 'default'
            ob.children = block.get('children', [])
            ob.frames = block.get('frames', [])
            obs_blocks1[ob_id] = ob
            # ignore frames for the moment

        obs_blocks2 = {}
        for key, obs in obs_blocks1.items():
            now = datetime.datetime.now()
            ob = ObservingBlock(instrument_id=obs.instrument, mode=obs.mode, start_time=now)
            ob.id = obs.uuid
            obs_blocks2[obs.id] = ob
            # FIXME: add alias, only if needed
            alias = ObservingBlockAlias(uuid=obs.uuid, alias=obs.id)

            # add frames
            # extract metadata from frames
            # FIXME:
            ingestdir = 'data'
            meta_frames = []
            for fname in obs.frames:
                full_fname = os.path.join(ingestdir, fname)
                print(fname, full_fname)
                result = metadata_fits(full_fname, self.drps)
                # numtype = result['type']
                # blck_uuid = obs.uuid # result.get('blckuuid', obs.uuid)
                result['path'] = fname
                meta_frames.append(result)

            for meta in meta_frames:
                # Insert into DB
                newframe = Frame()
                newframe.name = meta['path']
                newframe.uuid = meta['uuid']
                newframe.start_time = meta['observation_date']
                # No way of knowing when the readout ends...
                newframe.completion_time = newframe.start_time + datetime.timedelta(seconds=meta['darktime'])
                newframe.exposure_time = meta['exptime']
                newframe.object = meta['object']
                ob.frames.append(newframe)

            # set start/completion time from frames
            if ob.frames:
                ob.object = meta_frames[0]['object']
                ob.start_time = ob.frames[0].start_time
                ob.completion_time = ob.frames[-1].completion_time

            # Facts
            # add_ob_facts(session, ob, ingestdir)

            # raw frames insertion
            # for frame in frames:
            # call_event('on_ingest_raw_fits', session, frame, frames[frame])

            self.session.add(ob)
            self.session.add(alias)

        # processes children
        for key, obs in obs_blocks1.items():
            if obs.children:
                # get parent
                parent = obs_blocks2[key]
                for cid in obs.children:
                    # get children
                    child = obs_blocks2[cid]
                    parent.children.append(child)

        print('stage4')
        for key, obs in obs_blocks2.items():
            if obs.object is None:
                o1, s1, c1 = complete_recursive_first(obs)
                o2, s2, c2 = complete_recursive_last(obs)
                obs.object = o1
                obs.start_time = s1
                obs.completion_time = c2

        self.session.commit()


    def search_oblock_from_id(self, obsref):

        return search_oblock_from_id(self.session, obsref)

    def search_prod_obsid(self, ins, obsid, pipeline):
        """Returns the first coincidence..."""
        ins_prod = None # self.prod_table[ins]

        # search results of these OBs
        for prod in ins_prod.values():
            if prod['ob'] == obsid:
                # We have found the result, no more checks
                return StoredProduct(**prod)
        else:
            raise NoResultFound('result for ob %i not found' % obsid)

    def search_prod_req_tags(self, req, ins, tags, pipeline):
        return self.search_prod_type_tags(req.type, ins, tags, pipeline)

    def search_prod_type_tags(self, tipo, ins, tags, pipeline):
        """Returns the first coincidence..."""

        _logger.debug('query search_prod_type_tags type=%s instrument=%s tags=%s pipeline=%s',
                      tipo, ins, tags, pipeline)
        # drp = self.drps.query_by_name(ins)
        label = tipo.name()
        # print('search prod', tipo, ins, tags, pipeline)
        session = self.session
        # FIXME: and instrument == ins
        res = session.query(DataProduct).filter(DataProduct.type == label).order_by(DataProduct.priority.desc())
        _logger.debug('requested tags are %s', tags)
        for prod in res:
            pt = {}
            # TODO: facts should be a dictionary
            for key, val in prod.facts.items():
                pt[val.key] = val.value
            # print('product ', prod.id, 'tags', pt)
            _logger.debug('found value with id %d', prod.id)
            _logger.debug('product tags are %s', pt)

            if tags_are_valid(pt, tags):
                _logger.debug('tags are valid, return product, id=%s', prod.id)
                _logger.debug('content is %s', prod.content)
                # this is a valid product
                return StoredProduct(id=prod.id,
                                     content=load(tipo, os.path.join(self.basedir, prod.content)),
                                     tags=pt
                                     )
            _logger.debug('tags are in valid')
        else:
            _logger.debug('query search_prod_type_tags, no result found')
            msg = 'type %s compatible with tags %r not found' % (label, tags)
            raise NoResultFound(msg)

    def search_param_type_tags(self, name, tipo, instrument, mode, pipeline, tags):
        _logger.debug('query search_param_type_tags name=%s instrument=%s tags=%s pipeline=%s mode=%s', name, instrument, tags, pipeline, mode)
        session = self.session

        if isinstance(instrument, six.string_types):
            instrument_id = instrument
        else:
            instrument_id = instrument.name

        res = session.query(RecipeParameters).filter(
            RecipeParameters.instrument_id == instrument_id,
            RecipeParameters.pipeline == pipeline,
            RecipeParameters.name == name,
            RecipeParameters.mode == mode).one_or_none()
        _logger.debug('requested tags are %s', tags)
        if res is None:
            raise NoResultFound("No parameters for %s mode, pipeline %s", mode, pipeline)
        for value in res.values:
            pt = {}
            for f in value.facts.values():
                pt[f.key] = f.value
            _logger.debug('found value with id %d', value.id)
            _logger.debug('param tags are %s', pt)

            if tags_are_valid(pt, tags):
                _logger.debug('tags are valid, param, id=%s, end', value.id)
                _logger.debug('content is %s', value.content)
                # this is a valid product
                return StoredParameter(value.content)
        else:
            raise NoResultFound("No parameters for %s mode, pipeline %s", mode, pipeline)

    def assembly_instrument(self, keyval, date, by_key='name'):
        return assembly_instrument(self.components, keyval, date, by_key=by_key)

    def obsres_from_oblock(self, oblock, as_mode=None):
        # From dictdal

        from numina.core.oresult import ObservationResult

        # Internal copy
        obsres = oblock
        obsres.configuration = 'default'
        obsres.mode = as_mode or obsres.mode
        _logger.debug("obsres_from_oblock_2 id='%s', mode='%s' START", obsres.id, obsres.mode)

        try:
            this_drp = self.drps.query_by_name(obsres.instrument)
        except KeyError:
            raise ValueError('no DRP for instrument {}'.format(obsres.instrument))

        # Reserved names
        if obsres.mode in self._RESERVED_MODE_NAMES:
            selected_mode = None # null mode
        else:
            selected_mode = this_drp.modes[obsres.mode]

        if selected_mode:
            # This is used if we pass a option here or
            # if mode.build_ob_options is defined in the mode class
            # seems useless
            obsres = selected_mode.build_ob(obsres, self)
            # Not needed, all the information is obtained from
            # the requirements
            obsres = selected_mode.tag_ob(obsres)

        _logger.debug('assembly instrument model')
        key, date_obs, keyname = this_drp.select_profile(obsres)
        obsres.configuration = self.assembly_instrument(key, date_obs, keyname)
        obsres.profile = obsres.configuration

        auto_configure = True
        sample_frame = obsres.get_sample_frame()
        if auto_configure and sample_frame is not None:
            _logger.debug('configuring instrument model with image')
            obsres.profile.configure_with_image(sample_frame.open())
        else:
            _logger.debug('no configuring instrument model')
        return obsres

    def obsres_from_oblock_id(self, obsid, as_mode=None, configuration=None):
        # Search
        _logger.debug('query obsres_from_oblock_id with obsid=%s', obsid)
        oblock = self.search_oblock_from_id(obsid)
        print(oblock.instrument)
        return self.obsres_from_oblock(oblock, as_mode)

    def search_parameter(self, name, tipo, obsres, options=None):
        # returns StoredProduct
        instrument = obsres.instrument
        mode = obsres.mode
        tags = obsres.tags
        pipeline = obsres.pipeline

        if name in self.extra_data:
            value = self.extra_data[name]
            content = StoredParameter(value)
            return content
        else:
            return self.search_param_type_tags(name, tipo, instrument, mode, pipeline, tags)


    def search_product(self, name, tipo, obsres, options=None):
        # returns StoredProduct
        ins = obsres.instrument
        tags = obsres.tags
        pipeline = obsres.pipeline

        if name in self.extra_data:
            val = self.extra_data[name]
            content = load(tipo, val)
            return StoredProduct(id=0, tags={}, content=content)
        else:
            return self.search_prod_type_tags(tipo, ins, tags, pipeline)

    def search_result_relative(self, name, tipo, obsres, mode, field, node, options=None):
        # So, if node is children, I have to obtain
        session = self.session
        if node == 'children':
            print('obtain', field, 'from all the children of', obsres.taskid)
            res = session.query(DataProcessingTask).filter_by(id=obsres.taskid).one()
            result = []
            for child in res.children:
                # this can be done better...
                nodes = session.query(ReductionResult).filter_by(task_id=child.id).first()
                # this surely can be a mapping instead of a list

                for prod in nodes.values:
                    if prod.name == field:
                        st = StoredProduct(
                            id=prod.id,
                            content=load(DataFrameType(), os.path.join(self.basedir, prod.contents)),
                            tags={}
                        )
                        result.append(st)
                        break
            return result

        elif node == 'prev':
            print('obtain', field, 'from the previous node to', obsres.taskid)
            res = session.query(DataProcessingTask).filter_by(id=obsres.taskid).one()
            # inspect children of my parent
            parent = res.parent
            if parent:
                print([child for child in parent.children])
                raise NoResultFound
            else:
                # Im top level, no previous
                raise NoResultFound
        else:
            pass # print(dest, type, obsres, mode, field, node)
