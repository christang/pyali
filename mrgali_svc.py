import falcon
import simplejson
from mrgali import Alignment


class AlignmentService:
    def on_post(self, req, resp):
        body = req.stream.read()
        if not body:
            raise falcon.HTTPBadRequest('bad request', 'no sequences')

        try:
            doc = simplejson.loads(body)
        except KeyError:
            raise falcon.HTTPBadRequest('bad request', 'no sequences')
        except simplejson.JSONDecodeError:
            raise falcon.HTTPBadRequest('bad request', 'cannot parse')

        try:
            ref_ali = doc['ref']
            a = Alignment.from_reference(ref_ali)
        except KeyError:
            raise falcon.HTTPBadRequest('bad request', 'missing ref')

        try:
            seq_alis = doc['seqs']
            for i, seq_ali in seq_alis:
                a.merge(i, seq_ali)
        except KeyError:
            raise falcon.HTTPBadRequest('bad request', 'missing seqs')

        resp.body = simplejson.dumps({'result': str(a)})


api = falcon.API()
api.add_route('/merge', AlignmentService())
