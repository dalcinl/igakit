import numpy as np
try:
    from mayavi import mlab
    from mayavi.tools import helper_functions as _helper
    from mayavi.tools import sources as _sources
    from mayavi.tools import tools as _tools
    from tvtk.api import tvtk
except ImportError:
    from enthought.mayavi import mlab
    from enthought.mayavi.tools import helper_functions as _helper
    from enthought.mayavi.tools import sources as _sources
    from enthought.mayavi.tools import tools as _tools
    from enthought.tvtk.api import tvtk
try:
    from vtk.util import colors
except ImportError:
    class colors(object):
        black   = (0.00, 0.00, 0.00)
        grey    = (0.75, 0.75, 0.75)
        white   = (1.00, 1.00, 1.00)
        #
        red     = (1.00, 0.00, 0.00)
        green   = (0.00, 1.00, 0.00)
        blue    = (0.00, 0.00, 1.00)
        #
        yellow  = (1.00, 1.00, 0.00)
        cyan    = (0.00, 1.00, 1.00)
        magenta = (1.00, 0.00, 1.00)
        #
        orange  = (1.00, 0.50, 0.00)

show = mlab.show

figure = mlab.figure
close = mlab.close
save = mlab.savefig

title = mlab.title
xlabel = mlab.xlabel
ylabel = mlab.ylabel
zlabel = mlab.zlabel
colorbar = mlab.colorbar

points3d = mlab.points3d
plot3d   = mlab.plot3d
quiver3d = mlab.quiver3d

def _grid_data(LINES, SURFS, VOLS):
    point_list = []
    line_list  = []
    poly_list  = []
    for (x, y, z) in LINES+SURFS+VOLS:
        assert x.ndim == y.ndim == z.ndim
        assert x.shape == y.shape, "Arrays x and y must have same shape."
        assert y.shape == z.shape, "Arrays y and z must have same shape."
        #
        xyz = [x.ravel(), y.ravel(), z.ravel()]
        points = np.column_stack(xyz).ravel()
        points.shape = (-1, 3)
        lines = np.zeros((0, 2), dtype='l')
        polys = np.zeros((0, 4), dtype='l')
        #
        grid = np.arange(x.size, dtype='l').reshape(x.shape)
        if x.ndim == 1:
            p0 = grid[:-1].ravel()
            p1 = grid[+1:].ravel()
            verts = [p0,p1]
            lines = np.column_stack(verts).ravel()
            lines.shape = (-1, 2)
        elif x.ndim == 2:
            p0 = grid[:-1, :-1].ravel()
            p1 = grid[+1:, :-1].ravel()
            p2 = grid[+1:, +1:].ravel()
            p3 = grid[:-1, +1:].ravel()
            verts = [p0,p1,p2,p3]
            polys = np.column_stack(verts).ravel()
            polys.shape = (-1, 4)
        elif x.ndim == 3:
            p0 = grid[:-1, :-1, :-1].ravel()
            p1 = grid[+1:, :-1, :-1].ravel()
            p2 = grid[+1:, +1:, :-1].ravel()
            p3 = grid[:-1, +1:, :-1].ravel()
            p4 = grid[:-1, :-1, +1:].ravel()
            p5 = grid[+1:, :-1, +1:].ravel()
            p6 = grid[+1:, +1:, +1:].ravel()
            p7 = grid[:-1, +1:, +1:].ravel()
            verts = [p0,p1,p2,p3,
                     p4,p5,p6,p7,
                     p0,p1,p5,p4,
                     p2,p6,p7,p3,
                     p0,p3,p7,p4,
                     p1,p2,p6,p5]
            polys = np.column_stack(verts).ravel()
            polys.shape = (-1, 4)
        point_list.append(points)
        line_list.append(lines)
        poly_list.append(polys)
    offset = 0
    for points, lines, polys in zip(point_list, line_list, poly_list):
        lines += offset
        polys += offset
        offset += len(points)
    return (np.row_stack(point_list),
            np.row_stack(line_list),
            np.row_stack(poly_list))

class MGridSource(_sources.MlabSource):

    lines = _helper.List(_helper.Array, [])
    surfs = _helper.List(_helper.Array, [])
    vols  = _helper.List(_helper.Array, [])

    def reset(self, **traits):
        self.set(trait_change_notify=False, **traits)
        lines, surfs, vols = self.lines, self.surfs, self.vols
        points, lines, polys = _grid_data(lines, surfs, vols)
        if self.dataset is None:
            pd = tvtk.PolyData()
        else:
            pd = self.dataset
        pd.set(lines=None, polys=None)
        pd.set(points=points)
        pd.set(lines=lines, polys=polys)
        pd.point_data.scalars = points[:,-1].copy()
        pd.point_data.scalars.name = 'scalars'
        self.dataset = pd

    @classmethod
    def grid_source(MGridSource, lines=(), surfs=(), vols=(), **kwargs):
        convert = _sources.convert_to_arrays
        lines = [convert((x, y, z)) for (x, y, z) in lines]
        surfs = [convert((x, y, z)) for (x, y, z) in surfs]
        vols  = [convert((x, y, z)) for (x, y, z) in vols ]
        data_source = MGridSource()
        data_source.reset(lines=lines, surfs=surfs, vols=vols)
        name = kwargs.pop('name', 'GridSource')
        ds = _tools.add_dataset(data_source.dataset, name, **kwargs)
        data_source.m_data = ds
        return ds

class Grid3d(_helper.Mesh):
    lines = _helper.List(_helper.Array, [])
    surfs = _helper.List(_helper.Array, [])
    vols  = _helper.List(_helper.Array, [])
    _source_function = _helper.Callable(MGridSource.grid_source)

grid3d = _helper.document_pipeline(Grid3d())

del  MGridSource, Grid3d

# --


_resolution = { 1:256, 2:128 }
