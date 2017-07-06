import numpy as np
import gsd
import gsd.fl
import gsd.hoomd
import networkx as nx
from freud import box, density

class RDF(object):
    '''Class containing rdf functions'''
    def getRDFFromGSD(gsd_path,start_frame,end_frame,type1,type2,rmax=5.0,dr=0.01):
        rdf = density.RDF(rmax=rmax, dr=dr)
        rdf.resetRDF()
        f = gsd.fl.GSDFile(gsd_path, 'rb')
        t = gsd.hoomd.HOOMDTrajectory(f)
        n_frames = len(t)
        if start_frame>=0:
            print('RDF will be calulated from frame {} to {}. Available frames are from {} to {}.'.format(start_frame,end_frame,0,n_frames))
            for i in range(start_frame, end_frame):
                snapshot = t[i]
                sim_box = snapshot.configuration.box
                fbox = box.Box(Lx=sim_box[0], Ly=sim_box[1], Lz=sim_box[2])
                l_pos = snapshot.particles.position
                #print(len(l_pos),l_pos)
                # filter A particles.
                A_pos = l_pos[np.where(snapshot.particles.typeid == type1)]
                B_pos = l_pos[np.where(snapshot.particles.typeid == type2)]
                #print(len(A_pos),A_pos)
                # accumulate
                rdf.accumulate(fbox, A_pos, B_pos)
            # get the center of the histogram bins
            r = rdf.getR()
            # get the value of the histogram bins
            y = rdf.getRDF()
        else:
            print('WARNING: start frame cannot be less than 0')
        return r,y

