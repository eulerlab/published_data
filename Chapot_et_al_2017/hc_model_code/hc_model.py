import subprocess
import numpy as np
import pandas as pd
from io import StringIO
from io import BytesIO
import re
from PIL import Image
from PIL import ImageDraw
import networkx as nx

global __retsim_dir
global __vid_dir


def set_dir(retsim_directory,vid_directory=None):
    global __retsim_dir
    global __vid_dir
    __retsim_dir=retsim_directory
    if not vid_directory is None:
        __vid_dir=vid_directory

def run(retsim_dir=None,vid_dir=None,expt='hc_local_final',filename=None,errfilename=None,d=None,R=False,channel_default=True,vid_args=None,pov_args=None,pov_fn='temp',pov_in_fn=None,pov_image='temp.png',povres=3000,recolor=False,print_err=0,**kwargs):
    
    if retsim_dir is None:
        if '__retsim_dir' in globals():
            retsim_dir=__retsim_dir
        else:
            print('Retsim directory has to be specified')
            return
    if vid_dir is None:
        if '__vid_dir' in globals():
            vid_dir=__vid_dir

    ninfo      = kwargs.get('ninfo',None)
    nvalfile   = kwargs.get('nvalfile',None)
    rectype   = kwargs.get('rectype',None)
    stimtype   = kwargs.get('stimtype',None)
    node_scale = kwargs.get('node_scale',None)
    hb_nscale  = kwargs.get('hb_nscale',None)
    complam    = kwargs.get('complam',0.05)
    scone_id   = kwargs.get('scone_id',None)
    nrepeats   = kwargs.get('nrepeats',None)
    bg_inten   = kwargs.get('bg_inten',None)
    minten     = kwargs.get('minten',None)
    s_inten    = kwargs.get('s_inten',None)
    m_inten    = kwargs.get('m_inten',None)
    pulsedur   = kwargs.get('pulsedur',None)
    pulsegap   = kwargs.get('pulsegap',None)
    drm        = kwargs.get('drm',None)
    dcrm       = kwargs.get('dcrm',None)
    dri        = kwargs.get('dri',None)
    vnoise     = kwargs.get('vnoise',0)
    if channel_default:
        cadist     = kwargs.get('cadist',1e-3)
        caprox     = kwargs.get('caprox',3e-4)
        casoma     = kwargs.get('casoma',3e-4)
        kdist      = kwargs.get('kdist',1e-5)
        kprox      = kwargs.get('kprox',1e-5)
        ksoma      = kwargs.get('ksoma',1e-5)
        capdt      = kwargs.get('capdt',5e-6)
        cappx      = kwargs.get('cappx',5e-6)
        capsm      = kwargs.get('capsm',5e-6)
        capkdt     = kwargs.get('capkdt',5e-6)
        capkpx     = kwargs.get('capkpx',5e-6)
        capksm     = kwargs.get('capksm',5e-6)
    else:
        cadist     = kwargs.get('cadist',None)
        caprox     = kwargs.get('caprox',None)
        casoma     = kwargs.get('casoma',None)
        kdist      = kwargs.get('kdist',None)
        kprox      = kwargs.get('kprox',None)
        ksoma      = kwargs.get('ksoma',None)
        capdt      = kwargs.get('capdt',None)
        cappx      = kwargs.get('cappx',None)
        capsm      = kwargs.get('capsm',None)
        capkdt     = kwargs.get('capkdt',None)
        capkpx     = kwargs.get('capkpx',None)
        capksm     = kwargs.get('capksm',None)

    ipre        = kwargs.get('ipre',None)
    istart      = kwargs.get('istart',None)
    istop       = kwargs.get('istop',None)
    istep       = kwargs.get('istep',None)
    prestimdur  = kwargs.get('prestimdur',None)
    stimdur     = kwargs.get('stimdur',None)
    poststimdur = kwargs.get('poststimdur',None)
    stimtip     = kwargs.get('stimtip',None)
    
    tstep       = kwargs.get('tstep',None)
    setploti    = kwargs.get('setploti',None)


    


    arguments=['./retsim']

    if not expt: #print help if no experiment name is given
        p=subprocess.Popen(arg_list,cwd=retsim_dir,stderr=subprocess.PIPE,universal_newlines=True)
        output=p.communicate()
        print(output[1])
        return
        #return output[1] #maybe replace this by own help text

    arguments+=['--expt',expt]

    if ninfo!=None:       arguments+=['--ninfo',str(ninfo)]
    if nvalfile!=None:    arguments+=['--nvalfile',nvalfile]
    if rectype!=None:     arguments+=['--rectype',str(rectype)]
    if stimtype!=None:    arguments+=['--stimtype',str(stimtype)]
    if node_scale!=None:  arguments+=['--node_scale',node_scale]
    if hb_nscale!=None:   arguments+=['--hb_nscale',hb_nscale]
    if complam!=None:     arguments+=['--complam',str(complam)]
    if scone_id!=None:    arguments+=['--scone_id',str(scone_id)]
    if nrepeats!=None:    arguments+=['--nrepeats',str(nrepeats)]
    if bg_inten!=None:    arguments+=['--bg_inten',str(bg_inten)]
    if minten!=None:      arguments+=['--minten',str(minten)]
    if s_inten!=None:     arguments+=['--s_inten',str(s_inten)]
    if m_inten!=None:     arguments+=['--m_inten',str(m_inten)]
    if pulsedur!=None:    arguments+=['--pulsedur',str(pulsedur)]
    if pulsegap!=None:    arguments+=['--pulsegap',str(pulsedur)]
    if cadist!=None:      arguments+=['--cadist',str(cadist)]
    if caprox!=None:      arguments+=['--caprox',str(caprox)]
    if casoma!=None:      arguments+=['--casoma',str(casoma)]
    if kdist!=None:       arguments+=['--kdist',str(kdist)]
    if kprox!=None:       arguments+=['--kprox',str(kprox)]
    if ksoma!=None:       arguments+=['--ksoma',str(ksoma)]
    if capdt!=None:       arguments+=['--capdt',str(capdt)]
    if cappx!=None:       arguments+=['--cappx',str(cappx)]
    if capsm!=None:       arguments+=['--capsm',str(capsm)]
    if capkdt!=None:      arguments+=['--capkdt',str(capkdt)]
    if capkpx!=None:      arguments+=['--capkpx',str(capkpx)]
    if capksm!=None:      arguments+=['--capksm',str(capksm)]
    if drm!=None:         arguments+=['--drm',str(drm)]
    if dcrm!=None:        arguments+=['--dcrm',str(dcrm)]
    if dri!=None:         arguments+=['--dri',str(dri)]
    if vnoise!=None:      arguments+=['--vnoise',str(vnoise)]

    if ipre!=None:        arguments+=['--ipre',str(ipre)]
    if istart!=None:      arguments+=['--istart',str(istart)]
    if istop!=None:       arguments+=['--istop',str(istop)]
    if istep!=None:       arguments+=['--istep',str(istep)]
    if prestimdur!=None:  arguments+=['--prestimdur',str(prestimdur)]
    if stimdur!=None:     arguments+=['--stimdur',str(stimdur)]
    if poststimdur!=None: arguments+=['--poststimdur',str(poststimdur)]
    if stimtip!=None:     arguments+=['--stimtip',str(stimtip)]
    
    if tstep!=None:       arguments+=['--tstep',str(tstep)]
    if setploti!=None:    arguments+=['--setploti',str(setploti)]

    if R: #render povray image
        if not d:
            d=1
        arguments+=['-d',str(d),'-R','-1',pov_fn+'.pov']
        subprocess.call(arguments,cwd=retsim_dir)
        if recolor:
            with open(retsim_dir+pov_fn+'.pov','r') as f:
                pov_file=f.read()
            if stimtype==2:
                pov_file=pov_file.replace('nc.pov','nc_w.pov')
                pov_file=pov_file.replace('Cyan','Gray20')
                pov_file=pov_file.replace('Red','Gray20')
                pov_file=pov_file.replace('Blue','OrangeRed')
            else:
                pov_file=pov_file.replace('nc.pov','nc_w.pov')
                pov_file=pov_file.replace('Cyan','Gray30')
                pov_file=pov_file.replace('Red','Gray30')
                pov_file=pov_file.replace('Blue','Gray30')
                pov_file=pov_file.replace('Magenta','Silver')
            if scone_id==2:
                    pov_file=pov_file.replace('<-7.9355,0.517,4.49025>, 1.2 pigment {Green}','<-7.9355,0.517,4.49025>, 1.2 pigment {Blue}')
                    pov_file=pov_file.replace('<-7.9355,0.517,-0.50975>, 1.2 pigment {Green}','<-7.9355,0.517,-0.50975>, 1.2 pigment {Blue}')
                    pov_file=pov_file.replace('<-7.9355,0.517,-0.50975>,0.08 pigment {Green}','<-7.9355,0.517,-0.50975>,0.08 pigment {Blue}')
                    pov_file=pov_file.replace('<-7.9355,0.517,4.49025> pigment {Green}','<-7.9355,0.517,4.49025> pigment {Blue}')
            with open(retsim_dir+pov_fn+'.pov','w') as f:
                f.write(pov_file)

        pov_arg_list=['povray']
        if not pov_in_fn:
            pov_in_fn=pov_fn
        if pov_args:
            pov_arg_list+=pov_args+['+i'+pov_in_fn+'.pov','+o'+pov_image]
        else:
            pov_arg_list+=['+h1000','+w1000','+i'+pov_in_fn+'.pov','+o'+pov_image]
        subprocess.call(pov_arg_list,cwd=retsim_dir)
        im=Image.open(retsim_dir+'/'+pov_image)
        return im

    if d: #return image of the model
        arg_list+=['-d',str(d),'-v']
        if not vid_dir:
            vid_dir=retsim_dir.split('models')[0]+'bin/'
        p=subprocess.Popen(arguments,cwd=retsim_dir,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=False)
        output=p.communicate()
        err_out=output[1].decode('utf8')
        vid_arg_list=[vid_dir+'vid']
        if vid_args:
            vid_arg_list+=vid_args
        q=subprocess.Popen(vid_arg_list+['-c'],cwd=retsim_dir,stdin=subprocess.PIPE,stdout=subprocess.PIPE,universal_newlines=False)
        output2=q.communicate(input=output[0])
        im=Image.open(BytesIO(output2[0]))

        if print_err:
            return [im,err_out]
        return im
    elif filename: #write output to file
        arguments+=['-1',filename]
        if errfilename:
            arguments+=['-2',errfilename]
        subprocess.call(arguments,cwd=retsim_dir)
        return
    else: #standard case: return data as pandas dataframe
        p=subprocess.Popen(arguments,cwd=retsim_dir,stdout=subprocess.PIPE,stderr=subprocess.PIPE,\
                   universal_newlines=True)
        output=p.communicate()
        data=output[0]
        data=re.sub(' +',' ',data)
        data=re.sub(' \n','\n',data)
        try:
          data=pd.read_csv(StringIO(data),comment='#',delimiter=' ',header=None)
        except pd.io.common.EmptyDataError:
          print('No data output')
          data=[]
        if print_err==1:
            return [data,output[1]]
        return data


def drcable (im_size,x1, y1, z1, dia, x2, y2, z2, dia2, dscale=1, n1dia=0, n2dia=0, color=None, hide=False):
    """ Draw cable segments to show 2D morphologies, adapted from NeuronC"""
    
    im=Image.new('RGBA',tuple(im_size+1),(255,255,255,0))
    del_x=im_size[0]/2
    del_y=im_size[1]/2
    #    Draw a cable element in 2D projection of 3D rotation.  First, 
    # transform the element's 2 end points by a 3D rotation matrix, and
    # then draw the 2D projection of the element in a subframe, in a 
    # standard position (0,0) and orientation (horizontal).  The subframe is 
    # translated ("origin()") to one endpoint of the element,
    # and rotated ("rotate()") by an amount equal to the original
    # element's 3D orientation angle projected into 2D. 

    if(dia2<0):
        dia2 = dia

    if (dscale<0):
        dscale = -dscale;
        tdia = dia
        tdia2 = dia2
    else:
        tdia  = dia *  dscale
        tdia2 = dia2 * dscale
    r1 = n1dia * dscale / 2.0
    r2 = n2dia * dscale / 2.0
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    disto = np.sqrt(dx*dx + dy*dy + dz*dz)
    if (disto==0.0):
        disto = 1e-6;

    tx1=x1;ty1=y1;tz1=z1
    tx2=x2;ty2=y2;tz2=z2
    dy = ty2-ty1
    dx = tx2-tx1
    dist = np.sqrt(dx*dx + dy*dy)
    dratio = dist / disto # amount of foreshortening
    r1 *= dratio
    r2 *= dratio
    if (abs(dx)<1e-8): dx = .0001
    theta = np.arctan(dy/dx)
    if (dx<0): theta += np.pi

    if color==None: color='green'

    if (hide):
        return(im)

    else:
        draw=ImageDraw.Draw(im)
        if (r1+r2 < dist):       # if cable not inside sphere
            width  = tdia  / 2.0
            width2 = tdia2 / 2.0
            fill = 1
            if (width > width2): # if tapered
                x1 = r1+del_x
                y1 = -width2+del_y
                x2 = dist-r2+del_x
                y2 = -width+del_y
                draw.polygon([(x1,y1),(x2,y1),(x1,y2)],fill=color)
                y1 = width2+del_y
                y2 = width+del_y
                draw.polygon([(x1,y1),(x2,y1),(x1,y2)],fill=color)
            elif (width < width2):
                x1 = dist-r2+del_x
                y1 = -width+del_y
                x2 = r1+del_x
                y2 = -width2+del_y
                draw.polygon([(x1,y1),(x2,y1),(x1,y2)],fill=color)
                y1 = width+del_y
                y2 = width2+del_y
                draw.polygon([(x1,y1),(x2,y1),(x1,y2)],fill=color)

            if (width > width2): width = tdia2 / 2.0
            x1 = r1+del_x
            y1 = -width+del_y
            x2 = dist-r2+del_x
            y2 = width+del_y
            draw.rectangle([x1,y1,x2,y2],fill=color)
        del draw
        return(im.rotate(-theta/np.pi*180,translate=(tx1-del_x,ty1-del_y)))


def create_graph(morph):
    nodes = morph['node'].values
    pos = np.array([morph['x'].values, morph['y'].values, morph['z'].values]).T
    radius = morph['dia'].values
    region = morph['region'].values
    pid = morph['parent'].values
    
    # create node data (Kick whatever info you don't need)
    node_keys = ['pos', 'type', 'dia']
    node_data = list(zip(nodes, [dict(zip(node_keys, [pos[ix], region[ix], radius[ix]])) for ix in range(pos.shape[0])]))

    # create edge data
    n_ = nodes.tolist()
    parent_idx = [n_.index(pid[ix]) for ix in range(1,len(pid))]
    ec = np.sqrt(np.sum((pos[parent_idx] - pos[1:]) ** 2, axis=1)) # euclidean distance btw two connected nodes
    edge_keys = ['euclidean_dist']

    edge_data = list(zip(pid[1:], nodes[1:],
                         [dict(zip(edge_keys, [ec[ix]])) for ix in range(ec.shape[0])]))

    G = nx.Graph()
    G.add_nodes_from(node_data)
    G.add_edges_from(edge_data)
    
    return G


def tip_distances(morph,cone_tips):
    G = create_graph(morph)
    cone_distance = np.zeros((cone_tips.shape[0],cone_tips.shape[0]))
    for i in range(1,cone_distance.shape[0]):
        for j in range(0,i):
            cone_distance[i,j] = cone_distance[j,i] = nx.shortest_path_length(G,cone_tips[j],cone_tips[i], weight='euclidean_dist')
    return cone_distance


def nodes_to_tip(morph,tip):
    G = create_graph(morph)
    nodes_to_tip=nx.shortest_path(G,0,tip, weight='euclidean_dist')
    distance_to_tip=[]
    for node in nodes_to_tip:
        distance_to_tip.append(nx.shortest_path_length(G,node,tip,weight='euclidean_dist'))
    nodes_to_tip=np.array([nodes_to_tip,distance_to_tip]).T
    return nodes_to_tip
    
