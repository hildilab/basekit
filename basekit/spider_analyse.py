from __future__ import with_statement
from __future__ import division
import zipfile
import os
import json
import math
import itertools
import numpy as np
from utils import copy_dict
from utils.tool import _, _dir_init, PyTool, CmdTool, ScriptMixin
from utils.numpdb import NumPdb
from utils.mrc import get_mrc_header, getMrc
from utils.math import rmsd,tmscore
from pdb import *#PdbSplit
from dssp import *
from scipy import stats
from itertools import groupby
#from pdb import PdbEdit, NumpdbTest
from operator import itemgetter
DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "spider" )
SPIDER_CMD = "spider" 



            
# 2010 Cryo-EM Modeling Challenge: http://ncmi.bcm.edu/challenge


class MrcHeaderPrinter( PyTool ):
    args = [
        _( "mrc_file", type="file", ext="mrc" )
    ]
    def func( self, *args, **kwargs ):
        header = get_mrc_header( self.mrc_file )
        for name, value in zip(header._fields, header):
            print "%20s\t%s" % ( name, value )


class Spider( CmdTool, ScriptMixin ):
    args = [
        _( "script_file", type="file", ext="spi" )
    ]
    script_tmpl = None
    tmpl_dir = TMPL_DIR
    def _init( self, script_file, *args, **kwargs ):
        if script_file=="__tmpl__":
            self.script_file = self.outpath( self.script_tmpl )
        script_ext = "spi"
        data_ext = "cpv"
        # spider spi/cpv @box
        self.cmd = [
            SPIDER_CMD, 
            "%s/%s" % ( script_ext, data_ext ), 
            "@%s" % self.relpath( self.script_file, no_ext=True )
        ]


class SpiderShift( Spider ):
    """Shifts and converts mrc file and pdb and shifts the """
    args = [
        _( "mrc_file", type="file", ext="map" ),
        _( "pdb_file", type="file", ext="pdb" ),
    ]
    out = [
        _( "map_shift", file="shift.mrc" ),
        _( "edited_pdb_file", file="edited.pdb" )
    ]
    script_tmpl = "shift.spi"

    def _init( self, *args, **kwargs ):
        super(SpiderShift, self)._init( "__tmpl__" )

    def _pre_exec( self ):
        boxsizex=getMrc(self.mrc_file,'nx' )
        boxsizey=getMrc(self.mrc_file,'ny' )
        boxsizez=getMrc(self.mrc_file,'nz' )
        originx=getMrc(self.mrc_file,'nxstart')*-1
        originy=getMrc(self.mrc_file,'nystart')*-1
        originz=getMrc(self.mrc_file,'nzstart')*-1
        xorg=getMrc(self.mrc_file,'xorg')*-1
        yorg=getMrc(self.mrc_file,'yorg')*-1
        zorg=getMrc(self.mrc_file,'zorg')*-1
        print "nxstart", abs(getMrc(self.mrc_file,'nxstart')),abs(getMrc(self.mrc_file,'nystart')),abs(getMrc(self.mrc_file,'nzstart')),
        print "xorg", abs(getMrc(self.mrc_file,'xorg')),abs(getMrc(self.mrc_file,'yorg')),abs(getMrc(self.mrc_file,'zorg'))
        print "origin", originx,originy, originz
        size=getMrc(self.mrc_file,'xlen' )
        order=getMrc(self.mrc_file,'mapc' )
        pixelsize=(size/boxsizex)
       # print
        if originx!=0 and xorg==0 :
            shx = (originx -(boxsizex/2)) * pixelsize
            shy = (originy -(boxsizey/2)) * pixelsize
            shz = (originz -(boxsizez/2)) * pixelsize
        else:
            shx = ((xorg  / pixelsize)-(boxsizex/2))*pixelsize
            shy = ((yorg/ pixelsize)-(boxsizey/2))*pixelsize
            shz = ((zorg / pixelsize)-(boxsizez/2))*pixelsize
        self._make_script_file(    
            mrc_file=self.relpath( self.mrc_file ),
            order=order,
            pixelsize=pixelsize,
            boxsizex=boxsizex,
            boxsizey=boxsizey,
            boxsizez=boxsizez,
            originx=originx,
            originy=originy,
            originz=originz,
            shx  =shx,
            shy=shy,
            shz=shz
        )

        
        print shx
        print shy
        print shz
        PdbEdit( 
            self.pdb_file, shift= [shx, shy, shz]
        )


class SpiderConvert( Spider ):
    """Simple tool that converts mrc files to the spider format"""
    args = [
        _( "mrc_file", type="file", ext="mrc" )
    ]
    out = [
        _( "map_file", file="mapupload.cpv" )
    ]
    script_tmpl = "convert.spi"

    def _init( self, *args, **kwargs ):
        super(SpiderConvert, self)._init( "__tmpl__" )

    def _pre_exec( self ):
        #
        order=getMrc(self.mrc_file,'mapc' )
        print "order", order          
        self._make_script_file( 
            mrc_file=self.relpath( self.mrc_file ),order=order
        )


class SpiderPdbBox( PyTool ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "result_file", type="file" ),
        _( "mrc_file", type="file", ext="mrc" )
    ]
    out = [
        _( "edited_pdb_file", file="edited.pdb" )
    ]

    def func( self ):
        #mrc_dic={}
        #mrc_list=mrc_header( self.mrc_file)
        #for name, value in zip(mrc_list._fields, mrc_list):
         #   mrc_dic[name]=value
            #print mrc_dic[param]
        #param='nx'
        with open( self.result_file, "r") as rf:
            content = rf.readlines ()
            al = content[1].strip ()
            rs = al.split ()
            bs = float (rs [2])
            ol1 = float (rs [3])
            ol2 = float (rs [4])
            ol3 = float (rs [5])
            ps = float (rs [7])
            
            bbsx=getMrc(self.mrc_file, 'nx')#mrc_dic["nx"]
            bbsy=getMrc(self.mrc_file, 'ny')#mrc_dic["nx"]
            bbsz=getMrc(self.mrc_file, 'nz')
            obs = bs*ps
            x =  (ol1 - (bbsx/2)) * ps - 1 
            y =  (ol2 - (bbsy/2)) * ps - 1
            z =  (ol3 - (bbsz/2)) * ps - 1
            print [x,y,z, bs,bbsx]
        PdbEdit (self.pdb_file, box = [ x, y, z, obs, obs, obs ] )


class SpiderDeleteFilledDensities( Spider ):
    args = [
        _( "mrc_file", type="file", ext="mrc" ),
        _( "map_file", type="file", ext="cpv", help="boxil.cpv" ),
        _( "pdb_file", type="file", ext="pdb" ),
        _( "result_file", type="file", ext="cpv", help="ergebnisse.cpv"),
        _( "resolution", type="float", range=[1, 10],
            help="of the map_file" ),
        _( "res1", type="sele", help="resno:chain, i.e. 10:A" ),
        _( "res2", type="sele" )
    ]
    out = [
        _( "empty_map_file", file="usermap.cpv" )
    ]
    script_tmpl = "delete_filled_densities.spi"
    def _init( self, *args, **kwargs ):
        super(SpiderDeleteFilledDensities, self)._init( "__tmpl__" )
    def _pre_exec( self ):
        boxsizex=getMrc(self.mrc_file,'nx' )
        boxsizey=getMrc(self.mrc_file,'ny' )
        boxsizez=getMrc(self.mrc_file,'nz' )
        size=getMrc(self.mrc_file,'xlen' )
        pixelsize=(size/boxsizex)
        print self.res1
        LoopDelete(self.pdb_file,self.res1['chain'],self.res1['resno']+1,self.res2['resno']-1)
        self._make_script_file( 
            map_name=self.relpath( self.map_file, no_ext=True ), 
            pdb_file=self.relpath( 'noloop.pdb' ),
            result_file=self.relpath( self.result_file, no_ext=True ),
            pixelsize=pixelsize,
            resolution=self.resolution,
            boxsizex=boxsizex,
            boxsizey=boxsizey,
            boxsizez=boxsizez,
            tmp_dir=self.relpath( self.output_dir ) + os.sep
        )


class SpiderBox( Spider ):
    args = [
        _( "mrc_file", type="file", ext="mrc" ),
        _( "map_file", type="file", ext="cpv" ),
        _( "pdb_file", type="file", ext="pdb" ),
        _( "res1", type="sele", help="resno:chain, i.e. 10:A" ),
        _( "res2", type="sele" ),
        _( "length", type="int", range=[1, 30] ),
        #_( "pixelsize", type="slider", range=[1, 10], fixed=True ),
        _( "resolution", type="float", range=[1, 10],
            help="of the map_file" )
    ]
    out = [
        _( "var_file", file="variables.cpv" ),
        _( "box_file", file="ergebnisse.cpv" ),
        _( "box_map_file", file="boxil.cpv" )
    ]
    script_tmpl = "box.spi"
    def _init( self, *args, **kwargs ):
        super(SpiderBox, self)._init( "__tmpl__" )
    def _pre_exec( self ):
        boxsize=getMrc(self.mrc_file,'nx' )
        size=getMrc(self.mrc_file,'xlen' )
        pixelsize=(size/boxsize)
        coords1, coords2 = self._get_coords( 
            self.pdb_file, self.res1, self.res2 
        )
        self._make_variables_file(
            coords1, coords2, self.length, pixelsize, self.resolution
        )
        self._make_script_file( 
            map_name=self.relpath( self.map_file, no_ext=True ), 
            var_name=self.relpath( self.var_file, no_ext=True )
        )
    def _get_coords( self, pdb_file, res1, res2 ):
        npdb = NumPdb( pdb_file )
        return npdb.center( **res1 ), npdb.center( **res2 )
    def _make_variables_file( self, coords1, coords2, length, 
                              pixelsize, resolution ):
        variables = "1 9 %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %i %4.2f %4.2f" % (
            coords1[0], coords1[1], coords1[2],
            coords2[0], coords2[1], coords2[2],
            length, pixelsize, resolution
        )
        with open( self.var_file, "w" ) as fp:
            fp.write( variables )


class SpiderCropMap ( Spider ):
    args = [
   
    _( "map_file", type="file", ext="cpv" ),
    _( "pdb_file", type="file", ext="pdb" ),
    _( "pixelsize", type="float", range=[1, 10], fixed=True )
    
    ]
    out = [
        
    
    _( "box_file", file="ergebnisse.cpv" ),
    _( "box_map_file", file="boxil.cpv" )
    ]
    script_tmpl = "crop2.spi"
    def _init( self, *args, **kwargs ):
        super(SpiderCropMap, self)._init( "__tmpl__" )    
    def _pre_exec (self):
        npdb = NumPdb( self.pdb_file )
        psf=1/self.pixelsize
        schrumpf=npdb['xyz']*psf
        #te=
        ma=np.amax(schrumpf, axis=0)
        mi=np.amin(schrumpf,axis=0)
        middle=[(mi[0]+ma[0])/2,(mi[1]+ma[1])/2,(mi[2]+ma[2])/2]
        abstand=round(((mi[0]-ma[0])**2+(mi[1]-ma[1])**2+(mi[2]-ma[2])**2)**0.5,0)
        print ma,mi,middle,abstand
        self._make_script_file(
        x1=middle[0],
        y1=middle[1],
        z1=middle[2],
        map_name=self.relpath(self.map_file,no_ext=True),
        abstand=abstand,
        ps=self.pixelsize,
        )

class SpiderReConvert( Spider ):
    args = [
        _( "box_file", type="file", ext="cpv" ),
        _( "map_file", type="file", ext="cpv" ),
        _( "box_map_file", type="file", ext="cpv" ),
        _( "ori_map_file", type="file", ext="cpv" )
    ]
    out = [
        _( "mrc_file", file="reconvert.mrc" ),
        _( "mrc_ori_file", file="reconvertori.mrc" )
    ]
    script_tmpl = "recon.spi"
    def _init( self, *args, **kwargs ):
        super(SpiderReConvert, self)._init( "__tmpl__" )
    def _pre_exec( self ):
        self._make_script_file( 
            box_name=self.relpath( self.box_file, no_ext=True ),
            map_name=self.relpath( self.map_file, no_ext=True ),
            box_map_name=self.relpath( self.box_map_file, no_ext=True ),
            ori_map_name=self.relpath( self.ori_map_file, no_ext=True )
        )


class SpiderCrosscorrelation( Spider ):
    args = [
        _( "map_file", type="file", ext="cpv" ),
        _( "box_map_file", type="file", ext="cpv" ),
        _( "box_file", type="file", ext="cpv" ),
        _( "loop_file", type="file", ext="pdb" ),
        _( "linkerinfo", type="file", ext="txt" ),
        _( "max_loops", type="int", range=[0, 200], default=100 )
    ]
    out = [
        _( "loop_dir", dir="loops" ),
        _( "crosscorrel_file", file="crosscorrelation.cpv" ),
        _( "crosscorrel_json", file="crosscorrelation.json" )
    ]
    script_tmpl = "crosscorrelation.spi"
    def _init( self, *args, **kwargs ):
        if not self.max_loops:
            self.max_loops = False
        super(SpiderCrosscorrelation, self)._init( "__tmpl__" )
    def _pre_exec( self ):
        self._split_loop_file()
        self._make_script_file( 
            map_name=self.relpath( self.map_file, no_ext=True ), 
            box_map_name=self.relpath( self.box_map_file, no_ext=True ), 
            box_name=self.relpath( self.box_file, no_ext=True ), 
            loop_dir=self.relpath( self.loop_dir ) + os.sep,
            max_loops=self.max_loops or 999
        )
    def _post_exec( self ):
        self._make_crosscorrel_json( compact=True )
       
    def _split_loop_file( self ):
        PdbSplit( 
            self.loop_file,self.linkerinfo, output_dir=self.loop_dir, backbone_only=True, 
            max_models=self.max_loops, resno_ignore=[ 1000, 2000 ], zfill=3
        )
    def _make_crosscorrel_json( self, compact=False ):
        crosscorrel_dict = {}
        with open( self.crosscorrel_file, "r" ) as fp:
            for line in fp:
                d = line.split()
                if len(d)==3:
                    crosscorrel_dict[ int(d[0]) ] = float( d[2] )
        with open( self.crosscorrel_json, "w" ) as fp:
            if compact:
                json.dump( crosscorrel_dict, fp, separators=(',',':') )
            else:
                json.dump( crosscorrel_dict, fp, indent=4 )
                


class SpiderSidechainCorrelation ( Spider ) :
    args = [
        _( "map_file", type="file", ext="cpv" ),
        _( "pdb_file", type="file", ext="pdb" ),
        _( "pixelsize", type="float", range=[0, 10] ),
        _( "resolution", type="float", range=[0, 10] ),
        _( "residue", type="int", range=[0, 999] ),
        _( "chain", type="str" )
    ]
    out = [
        _( "sidechain_dir", dir="rotamere" ),
        _( "sccrosscorrel_file", file="sccrosscorrelation.cpv" ),
        _( "crosscorrel_json", file="crosscorrelation.json" )
    ]
    script_tmpl = "sidechaincc.spi"
    def _init( self, *args, **kwargs ):
        super(SpiderSidechainCorrelation, self)._init( "__tmpl__" )    
    def _pre_exec( self ):
        npdb=NumPdb( self.pdb_file )
        sele={"resno": self.residue, "chain": self.chain}
        resname1=npdb.get('resname', **sele)[0]
        num_rota=get_rotno(resname1)
        coords=self._get_ca(self.pdb_file,self.chain, self.residue)
        print "huh" ,coords [0] [0]

        x_ca= coords[0] [0]
        y_ca=coords [0][1]
        z_ca=coords [0] [2]
        self._make_script_file(
            resolution=self.resolution,
            pixelsize=self.pixelsize,
            map_name=self.relpath( self.map_file, no_ext=True ),
            x_ca=x_ca ,
            y_ca=y_ca,
            z_ca=z_ca,
            cc_dir=self.relpath( self.sidechain_dir ) + os.sep,
            num_rota=num_rota,
            resname1=self.relpath (self.sidechain_dir)+ os.sep+resname1,
            residue=self.residue
            )
        MakeAllRotameres(self.pdb_file, self.chain,self.residue,zfill=2,output_dir=self.sidechain_dir)
    def _get_ca (self,pdb_file,chain, residue):
        npdb=NumPdb( pdb_file )
        sele1={"resno": self.residue, "chain": self.chain}
        resname1=npdb.get('resname', **sele1)[0]
        sele={"resno": residue, "chain": chain, "resname": resname1 ,"atomname":'CA'}
        cacoord=npdb.get('xyz',**sele)
        return cacoord
    #def _make_rota ( self, pdb_file,chain,residue,zfill) :
    #    print "hallo"
    #    MakeAllRotameres(self.pdb_file, self.chain,self.residue,self.residue,zfill=2)
    #    
        #pass  
            
            
class SpiderDeleteBackbone (Spider):
    args = [
     _( "map_file", type="file", ext="cpv" ),
     _( "pdb_file", type="file", ext="pdb" ),
     _( "pixelsize", type="float", range=[0, 10] )
        
            
     ]
    out = [
        _( "edited_pdb_file", file="nobb.pdb" ),
        _( "delete_map_file", file="deletebb.cpv" )
                
    ]
    script_tmpl = "delete_backbone.spi"
    def _init( self, *args, **kwargs ):
        super(SpiderDeleteBackbone, self)._init( "__tmpl__" )    
    def _pre_exec( self ):
        backbone = ( ' N  ',' C  ', ' CA ',' O  ' )
        npdb=NumPdb( self.pdb_file )
        back=npdb.sele(atomname=backbone)
        npdb.copy(sele=back).write('nobb.pdb')
        self._make_script_file(
            map_name=self.relpath(self.map_file, no_ext=True) ,
            pixelsize=self.pixelsize,
            pdb_file='nobb.pdb'
        )


class LoopSidechainCorrelation (PyTool):            
    args = [
     _( "map_file", type="file", ext="cpv" ),
     _( "pdb_file", type="file", ext="pdb" ),
     _( "pixelsize", type="float", range=[0, 10] ),
     _( "resolution", type="float", range=[0, 10] ),    
            
     ]
    out = [
        _( "sidechain_dir", dir="sidechains" )           
    ]
    tmpl_dir = TMPL_DIR
    def _init( self, *args, **kwargs ):
        
        self.delete_backbone= SpiderDeleteBackbone(
            self.map_file,
            self.pdb_file,
            self.pixelsize,
            **copy_dict( kwargs, run=False )
        )
    def func( self ):
        self.delete_backbone()
                
        npdb=NumPdb( self.pdb_file )
        for i, numa in enumerate( npdb.iter_resno() ):
            sele={'resno':i+1}
            resname1=numa.get('resname') [0]
            resno1=numa.get('resno') [0]
            chain1=numa.get('chain') [0]
            dirg="%s_%s_%i" %(chain1,resname1,resno1)
            
            print resname1, resno1, chain1
            if   numa.get('resname') [0] not in ('ALA', 'GLY'):
                if resname1 in  ('SER','CYS','PRO','ASP','THR','ASN','VAL','GLU','GLN'):
                    SpiderSidechainCorrelation(self.map_file, self.pdb_file, self.pixelsize,self.resolution, resno1, chain1,output_dir=self.subdir(dirg))
                else:            
                    SpiderSidechainCorrelation('deletebb.cpv', self.pdb_file, self.pixelsize,self.resolution, resno1, chain1,output_dir=self.subdir(dirg))
                aaname="%s.%s" %(dirg,'pdb')
                print aaname
                aa=npdb.sele(resname=resname1,resno=resno1)
                dire="%s/%s" %(dirg,aaname)
                npdb.copy(sele=aa).write(dire)
                OriSidechainCorrel('deletebb.cpv',dire,self.pixelsize,self.resolution,resno1,chain1,output_dir=self.subdir(dirg))
            else:
                aaname="%s.%s" %(dirg,'pdb')
                print aaname
                aa=npdb.sele(resname=resname1,resno=resno1)
                dire="%s/%s" %(dirg,aaname)
                if not os.path.exists(self.subdir(dirg)): os.makedirs(self.subdir(dirg))
                npdb.copy(sele=aa).write(dire)

      
class OriSidechainCorrel ( Spider ):
    args = [
    _( "map_file", type="file", ext="cpv" ),  
     _( "pdb_file", type="file", ext="pdb" ),
     _( "pixelsize", type="slider", range=[0, 10], fixed=True ),
     _( "resolution", type="slider", range=[0, 10], fixed=True )  ,
     _( "residue", type="slider", range=[0, 999],fixed=True),
    _( "chain", type="str")
     ]
    out = [
        _( "sidechain_dir", dir="sidechains" )           
    ]
    script_tmpl = "orisidechain.spi"
    def _init( self, *args, **kwargs ):
        super(OriSidechainCorrel, self)._init( "__tmpl__" )    
    def _pre_exec( self ):
        npdb=NumPdb( self.pdb_file )
        sele={"resno": self.residue, "chain": self.chain}
        resname1=npdb.get('resname', **sele)[0]
        num_rota=get_rotno(resname1)
        coords=self._get_ca(self.pdb_file,self.chain, self.residue)
        print "huh" ,coords [0] [0]
        x_ca= coords[0] [0]
        y_ca=coords [0][1]
        z_ca=coords [0] [2]
        self._make_script_file(
        x_ca= x_ca,
        y_ca=y_ca,
        z_ca=z_ca,
        map_name=self.relpath(self.map_file, no_ext=True) ,
        pixelsize=self.pixelsize,
        resolution=self.resolution,
        pdb_file=self.relpath(self.pdb_file, no_ext=False) ,
        )
    #def _post_exec( self ):
    #    self._make_result_pdb()

    def _get_ca (self,pdb_file,chain, residue):
        npdb=NumPdb( pdb_file )
        sele1={"resno": self.residue, "chain": self.chain}
        resname1=npdb.get('resname', **sele1)[0]
        sele={"resno": residue, "chain": chain, "resname": resname1 ,"atomname":'CA'}
        cacoord=npdb.get('xyz',**sele)
        return cacoord
##build the pdb file from the best rotamers


class OptimizeRotamer ( PyTool ):
    args = [
    _("result_direc",type="str")
    ]
    out=[
    _( "verybestrotamers", file="verybestrotamers.pdb" )
    
    ]
        
    def  _init( self , *args, **kwargs ):
        self.bestbuild=BuildBest(self.result_direc)
    def func (self):
        self.bestbuild()
        npdb=NumPdb( "bestrotamers.pdb" )
        
        a=True
        z=2
        #while a==True:
        clashes,a=find_all_clashes(npdb)
        gtree=get_tree(npdb['xyz'])
        print clashes
        for i in clashes:
            localresi=[]
            #hole CA atom vom ersten clashpartner 
            r=npdb.sele(resno=i[0],atomname='CA')
            p=npdb.get('xyz',sele=r)
            #indices im 6 A umkreis um den clash
            l=gtree.query_ball_point(p,6)
            print l[0]
            #aus den indices eine liste mit residuname und nummer machen
            for x in l[0]:
                a=npdb.get('resname')[x]
                b=npdb.get('resno')[x]
                c=npdb.get('chain')[x]
                localresi.append([b,a,c])
            localresi=sorted(localresi)
            localresi=list(localresi for localresi,_ in itertools.groupby(localresi))
            #print localresi
            for index,i in enumerate (localresi):
                if i[1] not in ('ALA','GLY'):
                    rotadir="%s_%s_%i" % (i[2],i[1],i[0])
                    ccsortpath=os.path.join(self.result_direc,rotadir)
                    ccsort="%s/%s" % (ccsortpath,'ccsort.cpv')
                    print ccsortpath
                    
            break

#class SpiderMinMap ( Spider ):
#    args = [
#    _("pdb_file", type="file"),  
    
class  LoopRotamerOptimize ( PyTool ):
    args = [
    _("pdb_file", type="file"),  
    _("result_direc",type="str"),
    _("residue_1",type="int"),
    _("residue_2",type="int")
    ]
    out=[
    _( "verybestrotamers", file="verybestrotamers.pdb" )
    
    ]
   # clashes fuer das ganze Protein bestimmen
    def func (self):
        npdb=NumPdb( self.pdb_file )
        clashes,a=find_all_clashes(npdb)
        gtree=get_tree(npdb['xyz'])
        lclashes=[]
        print clashes
        resi=range(self.residue_1,self.residue_2)
        print resi
        for i in clashes:
            print i[0], i[1]
            if (i[0] or i[1]) in resi:
                lclashes.append(i)
            else:
                continue
        for x in lclashes:
            print x 
        
            as1=x[0]
            as2=x[1]
            resname100=npdb.get('resname',resno=as1)[0]
            resname101=npdb.get('resname',resno=as2)[0]
            #    #print resname100, resname101, i
            if (resname100 and resname101) not in ('ALA','GLY'):
                    chain100=npdb.get('chain',resno=as1)[0]
                    chain101=npdb.get('chain',resno=as2)[0]
                    #print resname101
                    resno1=get_rotno(resname100)
                    resno2=get_rotno(resname101)
                    if resno1 >= resno2 :
            #            #print "anzahl rota" ,resno1,resno2
                        sele={'resno':as1}
                        atom_no= npdb.index(**sele)
            #            #erster clashpartner
                        rotadir="%s_%s_%i" % (chain100,resname100,as1)
                        ccsortpath=os.path.join(self.result_direc,rotadir)
                        ccsort="%s/%s" % (ccsortpath,'ccsort.cpv')
                    else:
                        sele={'resno':as2}
                        atom_no= npdb.index(**sele)
                        rotadir="%s_%s_%i" % (chain101,resname101,as2)
                        ccsortpath=os.path.join(self.result_direc,rotadir)
                        ccsort="%s/%s" % (ccsortpath,'ccsort.cpv')
            #        
            #        with open (ccsort,"r") as fil:
            #            file_lines=fil.readlines()
            #            loc_file = file_lines[z].split()
            #            nbr=str(loc_file [0]).zfill(2)
            #            rota="%s/%s/%s_%s.%s" % (ccsortpath,'rotamere',rotadir[2:],nbr,'pdb')
            #        newrot=NumPdb(rota)
            #        newroti=newrot.index()
            #        all_coords = npdb['xyz']
            #        for index, at_num in enumerate(atom_no):
            #            all_coords[at_num] = newrot['xyz'][index]
            #        npdb['xyz'] = all_coords
            #    else:
            #        print "buhu"
            #print "round", z-1
            #z +=1
            #npdb.write("bestrotamers.pdb")
            #roundd="%s_%i.%s" % ("bestrotaround",z-1,'pdb')
            #npdb.write(roundd)
            #clashes,a=find_all_clashes(npdb)

        
class BuildBest ( PyTool ):
    args = [
        _("result_direc",type="str")
        ]
    out=[
        _( "bestrotamers", file="bestrotamers.pdb" )
        
    ]
    def func (self):
        b=open(self.bestrotamers, "w")
        gesamt=[]
        for fn in os.listdir(self.result_direc):
            #print fn[2:5]
            path=os.path.join(self.result_direc,fn)
            ccsort="%s/%s" % (path, 'ccsort.cpv')
            if fn[2:5] not in ('ALA','GLY'):
                with open (ccsort,"r") as fil:
                    file_lines=fil.readlines()
                    loc_file = file_lines[1].split()
                    nbr=str(loc_file [0]).zfill(2)
                    rota="%s/%s/%s_%s.%s" % (path,'rotamere',fn[2:],nbr,'pdb')
                with open (rota, "r") as pd:
                    pdb_lines=pd.readlines()
                    
                    for lines in pdb_lines:
                        if lines.startswith ("END"):
                            continue
                        else:
                            gesamt.append(lines)
            else:
                norota="%s/%s.%s" % (path,fn,'pdb')
                with open (norota, "r") as pd:
                    pdb_lines=pd.readlines()
                    
                    for lines in pdb_lines:
                        if lines.startswith ("END"):
                            continue
                        else:
                            gesamt.append(lines)   
        gesamt.sort(key=lambda x: x[7:11])              
        
        b.writelines(gesamt)


class SideChainStatistics ( PyTool ):
        args = [
    _( "dataset_dir", type="dir" )
    ]
        out = [
     _( "stat_file", file="statistics.csv" ),
     _( "res_file", file="overall.txt" )
    ]
        #def rmsd( coords1, coords2 ):
        #    return np.sqrt( 
        #        np.sum( np.power( coords1-coords2, 2 ) ) / coords1.shape[0]
        #    )  
        def func (self):
            
            with open (self.stat_file, 'w') as rf:
                title="%s,%s,%s,%s,%s,%s,%s%s" %  ("residue","rotanumber","hmm","corr", "oricorr","rmsd","bigger","\n")
                rf.write(title)
                fi=sorted((os.listdir(self.dataset_dir)), key= lambda x: int(x.split('_')[-1]))
                g=0
                b=0
                m=0
                a=0#fi.sort(key=float)
                #print fi
                for pn in fi:#os.listdir("chainA"):
                   # print pn
                    #if os.path.isdir(pn):
                    z=0
                    #print pn
                    if pn[2:5]  not in ("GLY", "ALA"):
                        a+=1
                        ccsort=  os.path.join(self.dataset_dir,pn,"ccsort.cpv")
                        ori=os.path.join(self.dataset_dir,pn,"oricrosscorrelation.cpv")
                        
                        #print ori
                        
                        with open (ccsort, 'r') as hu:
                            lines=hu.readlines()
                            hurz=','.join(lines[1].split())
                            #hurz2=hurz+"\n"
                            ncor=lines[1].split()[-1]
                            
                        with open (ori,'r') as bu:
                            line2=bu.readlines()
                            #hurz.append(line2)
        
                            corr=line2[1].split()[-1]
                            #print ncor
                            if float(ncor)>float(corr):
                                z=1
                            
                        oripdb=os.path.join(self.dataset_dir,pn,pn+".pdb")
                        
                        best=(lines[1].split()[0]).zfill(2)
                        #print best
                        bestrota="%s_%s.%s" % (pn[2:],best,'pdb')
                        #print bestrota
                        rotadir=os.path.join(self.dataset_dir,pn,"rotamere")
                        bestpdb=os.path.join(self.dataset_dir,pn,"rotamere",bestrota)

                        npdb = numpdb.NumPdb( bestpdb, {
                            "phi_psi": False,
                            "sstruc": False,
                            "backbone_only": False,
                            "protein_only": False,
                            "detect_incomplete": False,
                            "configuration": False,
                            "info": False
                            })
                        sele = npdb.sele()
                        #npdb2=get_npdb( oripdb )
                        npdb2 = numpdb.NumPdb( oripdb, {
                            "phi_psi": False,
                            "sstruc": False,
                            "backbone_only": False,
                            "protein_only": False,
                            "detect_incomplete": False,
                            "configuration": False,
                            "info": False
                            })
                        sele2 = npdb2.sele()
                        a1=np.array(npdb['xyz'])
                        a2=np.array(npdb2['xyz'])
                        #print npdb['xyz']
                        rootm=rmsd(npdb['xyz'],npdb2['xyz'])
                        bestrmsd=[]
                        for file in os.listdir(rotadir):
                            
                            if file.endswith(".pdb"):
                                #print file, "hallo"
                                npdb3 = numpdb.NumPdb(os.path.join(rotadir,file), {
                            "phi_psi": False,
                            "sstruc": False,
                            "backbone_only": False,
                            "protein_only": False,
                            "detect_incomplete": False,
                            "configuration": False,
                            "info": False
                            })
                                #print npdb3['xyz']
                                rootm2=rmsd(npdb3['xyz'],npdb2['xyz'])
                                bestrmsd.append(float(rootm2))
                        test100=sorted(bestrmsd)
                        print file, test100 [0], test100 [-1]
                        if test100[0]==rootm:
                            print "yeahaa!"
                            y=1
                            g+=1
                        if rootm <1.2 and not test100[0]==rootm:
                            m+=1
                            y=2
                            
                        if rootm >=1.2 and not test100[0]==rootm:
                            y=0
                            b+=1
                            
                        hurz2="%s,%s,%s,%s,%s,%s%s" %  (pn,hurz,corr,rootm,z,y,"\n")
                        rf.write(hurz2)
                with open (self.res_file, 'w') as gh:
                    juhu="%s %s %s %s" % (a,g,m,b)
                    gh.write(juhu)
                print a,g,m,b
     
class SpiderCropMrc ( PyTool ):
    args = [
   
    _( "map_file", type="file", ext="mrc" ),
    _( "pdb_file", type="file", ext="pdb" ),
    _( "pixelsize", type="float", range=[1, 10], fixed=True )
    
    ]
    out = [
        _( "cropped_mrc", file="cropped.mrc" )
    ]
    tmpl_dir = TMPL_DIR
    def _init( self, *args, **kwargs ):
        self.spider_convert= SpiderConvert(
                                self.map_file,
                                **copy_dict(
                                    kwargs, run=False, output_dir=self.subdir("convert")
                                )
        )
        self.spider_crop= SpiderCropMap(
                            self.spider_convert.map_file,self.pdb_file,self.pixelsize,
                            **copy_dict(
                                    kwargs, run=False, output_dir=self.subdir("crop")
                                )
        )
        self.spider_reconvert=SpiderReConvert(
                                    self.spider_crop.box_file,
                                    self.spider_convert.map_file,
                                    self.spider_crop.box_map_file,
                                    self.spider_convert.map_file,
                                    **copy_dict(
                                    kwargs, run=False, output_dir=self.subdir("mrc")
                                )
                                    
        )
    def func( self ):
        self.spider_convert()
        self.spider_crop()
        self.spider_reconvert()
class Calcoricc( Spider ):
    args = [
    _( "boxsize",type="float"),
    _( "boxsize2",type="float"),
    _( "boxsize3",type="float"),
    _( "oriboxsize",type="float"),
    _( "map_file", type="file", ext="cpv" ),
    _( "mask_file",type="file", ext="cpv" ),
    _( "pdb_file", type="file", ext="pdb" ),
    _( "ecke1", type="float"),
    _( "ecke2", type="float"),
    _( "ecke4" , type="float"),
    _( "resolution", type="float"),
    _( "pixelsize", type="float")

    
    ]
    out = [
     _( "stat_file", file="statistics.csv" )
     
    ]
    #tmpl_dir = TMPL_DIR
    script_tmpl = "oricc.spi"
    def _init ( self, *args, **kwargs ):
        super(Calcoricc, self)._init( "__tmpl__" )   
    def _pre_exec( self ):
        self._make_script_file(
            map_file=self.relpath(self.map_file,no_ext=True),
            mask_file=self.relpath(self.mask_file,no_ext=True),
            pdb_file=self.relpath(self.pdb_file),
            ecke1=self.ecke1,
            ecke2=self.ecke2,
            ecke4=self.ecke4,
            bs1=self.boxsize,
            bs2=self.boxsize2,
            bs3=self.boxsize3,
            obs=self.oriboxsize,
            ps=self.pixelsize,
            rs=self.resolution
            )
class BuildModel( PyTool ):
    args = [
        _( "fragment_length", type="int"),
        _( "dataset_dir", type="dir" ),
        _( "oripdb", type="file")
    ]
    out = [
       _( "pdb_file", file="full_modell.pdb" )
    ]
    def func ( self, *args, **kwargs ):
        with open (self.pdb_file, 'w') as pdbmodel:
            
            chainfolders=os.listdir(self.dataset_dir)
            for g in chainfolders:
                fragfolders=os.path.join(self.dataset_dir,g,str(self.fragment_length))
                if os.path.isdir(fragfolders):
                    fi=sorted((os.listdir(fragfolders)), key= lambda x: int(x.split('_')[-1]))
                    hu=fi#[::1]
                    #print len(fi), len(hu)
                    #print hu
                    start=int(hu[0].split('_')[0])
                    start+=self.fragment_length
                    print start
                    for i in hu:
                        print i, start
                        ccsort=os.path.join(fragfolders,i,'loop_correl/crosscorrelation/ccsort.cpv')
                        if os.path.isfile(ccsort) and int(i.split('_')[0])==start:
                            print i
                            with open (ccsort,'r') as hurz:
                                lines=hurz.readlines()
    
                                nr=lines[1].split()[0]
                                corr=lines[1].split()[2]
                                name="%s_%s.%s" % (nr.zfill(3),'bb','pdb')
    
                                bestloop=os.path.join(self.dataset_dir,g,str(self.fragment_length),i,'loop_correl','crosscorrelation','loops',name)
                                #print bestloop
                                with open (bestloop, 'r') as bloopfile:
                                    for i, line in enumerate( bloopfile ):
                                        if line.startswith("ATOM"):
                                            pdbmodel.write(line)
                            start+=self.fragment_length
                        if not os.path.isfile(ccsort) and int(i.split('_')[0])==start:
                        #else:
                            start+=self.fragment_length
                        #    #print start
            hf="END"
            pdbmodel.write(hf)
            #print self.oripdb
            npdb=numpdb.NumPdb(self.oripdb)
            fmodel=numpdb.NumPdb( 'full_modell.pdb' )
            #rootm2=rmsd(npdb['xyz'],fmodel['xyz'])
            #print 'rmsd', rootm2
            
class Looporicc( PyTool ):
    args = [
        _( "boxsize",type="float"),
        _( "boxsize2",type="float"),
        _( "boxsize3",type="float"),
        _("pdb_file", type="file"),
        _( "dataset_dir", type="dir" )
        ]
    def func ( self, *args, **kwargs ):
            #l=25
        for l in range(5,36,1):
            #l=19
            chainfolders=os.listdir(self.dataset_dir)
            for g in chainfolders:
                
                loopfolder=os.path.join(self.dataset_dir,g,str(l))
                
                #os.path.exists(file_path)
                if os.path.isdir(loopfolder):
                    print 'G',g
                    fi=sorted((os.listdir(loopfolder)), key= lambda x: int(x.split('_')[-1]))#code
                    for x in fi:
                        name="%s.%s" % (x,'pdb')
                        oriccfile=os.path.join(loopfolder,x,'oricrosscorrelation.cpv')
                        if os.path.exists(oriccfile):
                            continue
                        else:
                            res1=x.split('_')[0]
                            res2=x.split('_')[1]
                            print res1,res2
                            npdb=numpdb.NumPdb( self.pdb_file)
                            test=npdb.copy(resno=[int(res1),int(res2)],chain=g)
                            out=os.path.join(loopfolder,x)
                            namefull= "%s/%s_%s.%s" % (out,res1,res2,'pdb')
                            test.write(namefull)
                            ori=os.path.join(loopfolder,x,name)
                            mapori=os.path.join(loopfolder,x,'loop_correl','delete_filled_densities','usermap.cpv')
                            mask=os.path.join(loopfolder,x,'loop_correl','crosscorrelation','mask.cpv')
                            #print maske
                            #print ori
                            shiftfile=os.path.join(loopfolder,x,'loop_correl','shift','shift.cpv')
                            eckenfile=os.path.join(loopfolder,x,'loop_correl','box','ergebnisse.cpv')
                            b_out=os.path.join(loopfolder,x,name)
                            
                            outd=os.path.join(loopfolder,x,'loop_correl')
                            ccsortf=os.path.join(loopfolder,x,'loop_correl','crosscorrelation','ccsort.cpv')
                            pdb_sh=os.path.join(loopfolder,x,'loop_correl','edited.pdb')
                            if os.path.isfile(eckenfile):
                               # print eckenfile
                                if os.path.isfile(ccsortf):
                                    print 'calculate',loopfolder,x,name
                                    with open (eckenfile, 'r') as ecki,open (shiftfile,'r') as shifti:
                                        guja=ecki.readlines()
                                        h=guja[1].split()
                                        obenlinks1=h[3]
                                        obenlinks2=h[4]
                                        obenlinks4=h[6]
                                        ps=h[7]
                                        reso=h[8]
                                        obs=h[2]
                                        #print h[4]
                                       # print out
                                        #print
                                        blue=shifti.readlines()
                                        f=blue[1].split()
                                        print f[2]
                                        sh1=float(f[2])
                                        sh2=float(f[3])
                                        sh3=float(f[4])
                                        PdbEdit(ori,shift=[sh1,sh2,sh3],output_dir=outd)
                                        
                                    Calcoricc(self.boxsize,self.boxsize2,self.boxsize3,obs,mapori,mask,pdb_sh,obenlinks1,obenlinks2,obenlinks4,reso,ps,output_dir=out)
                                    os.remove(b_out)
                                    os.rename(pdb_sh,b_out)    
class LoopCrosscorrel( PyTool ):
    args = [
        _( "mrc_file", type="file", ext="mrc" ),
        _( "pdb_file", type="file", ext="pdb" ),
        _( "loop_file", type="file", ext="pdb" ),
        _( "linkerinfo", type="file", ext="txt" ),
        _( "res1", type="sele" ),
        _( "res2", type="sele" ),
        _( "length", type="int", range=[1, 30] ),
        _( "resolution", type="float", range=[1, 10], step=0.1 ),
        _( "max_loops", type="int", range=[0, 500], default=100 )
    ]
    out = [
        _( "cropped_pdb", file="cropped.pdb" ),
        _( "ori_pdb_linker_file3", file="ori_pdb_linker_file3.pdb")
    ]
    tmpl_dir = TMPL_DIR
    def _init( self, *args, **kwargs ):
        #self.spider_shift = SpiderShift(
        #    self.mrc_file,self.pdb_file,
        #    **copy_dict( 
        #        kwargs, run=False, output_dir=self.subdir("shift") 
        #    )
        #)
        #self.spider_convert = SpiderConvert( 
        #    self.spider_shift.map_shift, 
        #    **copy_dict( 
        #        kwargs, run=False, output_dir=self.subdir("convert") 
        #    )
       # )
        self.spider_box = SpiderBox(
            self.mrc_file,
            #self.spider_shift.map_shift,
            "/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/1/files/mapupload.cpv",
            #self.spider_convert.map_file, 
            #self.spider_shift.edited_pdb_file
            "/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/1/files/edited.pdb", self.res1, self.res2, 
            self.length, self.resolution,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("box") )
        )
        self.pdb_box =  SpiderPdbBox(
           #self.spider_shift.edited_pdb_file,
           "/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/1/files/edited.pdb",
           self.spider_box.box_file,
           self.mrc_file,
           #self.spider_shift.map_shift,
           **copy_dict( kwargs, run=False, output_dir=self.subdir("pdbbox"))
        )
        self.spider_delete_filled_densities = SpiderDeleteFilledDensities(
            self.mrc_file,self.spider_box.box_map_file,self.pdb_box.edited_pdb_file,self.spider_box.box_file,
            
            #self.spider_shift.map_shift,self.spider_box.box_map_file, self.pdb_box.edited_pdb_file, self.spider_box.box_file,
            self.resolution,self.res1,self.res2,
            **copy_dict( 
                kwargs, run=False, 
                output_dir=self.subdir("delete_filled_densities") 
            )
        )
        self.spider_reconvert = SpiderReConvert(
            self.spider_box.box_file,
            "/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/1/files/mapupload.cpv",
            #self.spider_convert.map_file,
            self.spider_box.box_map_file,
            "/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/1/files/mapupload.cpv",
            #self.spider_convert.map_file,
            **copy_dict( 
                kwargs, run=False, output_dir=self.subdir("reconvert") 
            )
        )
        self.spider_crosscorrelation = SpiderCrosscorrelation(
            "/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/1/files/mapupload.cpv",
            #self.spider_convert.map_file, 
            self.spider_delete_filled_densities.empty_map_file, 
            self.spider_box.box_file, 
            self.loop_file,
            self.linkerinfo,
            **copy_dict( 
                kwargs, run=False, output_dir=self.subdir("crosscorrelation"),
                max_loops=self.max_loops
            )
        )
        self.output_files.extend( list( itertools.chain(
            #self.spider_convert.output_files, 
            self.spider_delete_filled_densities.output_files,
            self.spider_box.output_files,
            self.spider_reconvert.output_files,
            self.spider_crosscorrelation.output_files
        )))
    def func( self ):
        self._crop_pdb()
        #self.spider_shift()
        #self.spider_convert()
        self.spider_box()
        self.pdb_box ()
        self.spider_delete_filled_densities()
        self.spider_reconvert()
        self.spider_crosscorrelation()
    
    def _post_exec( self ):
        self.backshift_linker()
        
    def _crop_pdb( self ):
        npdb = NumPdb( self.pdb_file )
        npdb.write( 
            self.cropped_pdb, 
            chain=self.res1["chain"], 
            resno=[ self.res1["resno"]+1, self.res2["resno"]-1 ],
            invert=True
        )
    
    def backshift_linker ( self ) :
        #shiftpath=self.relpath(self.spider_shift.output_dir)
        shiftfile='/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/1/files/shift.cpv'
        #shiftfile="%s/%s" % (shiftpath,'shift.cpv')
        with open (shiftfile, "r") as rs:
            file_lines=rs.readlines()
            loc_file = file_lines[1].strip()
            shxb=float(loc_file[5:13])*-1
            shyb=float(loc_file[19:27])*-1
            shzb=float(loc_file[33:42])*-1

            backshift= np.array([shxb,shyb,shzb])
   
        loop_dir=self.relpath(self.spider_crosscorrelation.loop_dir)
        outputdir=self.subdir("oriloops")
        with open ( self.linkerinfo , 'r') as li:
            lilines=li.readlines()
        with open (self.loop_file, 'r') as lf:
            with open (self.ori_pdb_linker_file3, 'w') as fp_out:
                for line in lf:
                    if line.startswith("ATOM"):
                        x = "%4.3f" % (float(line[30:37]) + shxb)
                        a = "%8s" % x

                        y = "%4.3f" % (float(line[38:45]) + shyb)
                        b = "%8s" % y
                        z = "%4.3f" % (float(line[46:53]) + shzb)
                        c = "%8s" % z
                        print b
                        line = line = line[0:29] + a + line[38:]
                        line = line = line[0:38] + b + line[46:]
                        line = line = line[0:46] + c + line[54:]

        
                        fp_out.write( line )
                    else:
                        continue
                fp_out.write( "END" )
        for i in os.listdir(loop_dir):
            if i.endswith(".pdb"):
                loopfile="%s/%s" % (loop_dir,i)
                
                outfile="%s/%s" %(outputdir,i)
                with open (loopfile, 'r') as lp:
                    lonr=int(i[0:3])
                    lino=(lonr*4)-2

                    title= "%s     %s %s ,%s ,%s ,%s ,%s" % ("TITLE","linker",lonr,lilines[lino].strip(),lilines[lino+1].strip(),lilines[lino+2].strip(),lilines[lino+3])

                    with open (outfile, 'w') as of:
                        of.write(title)
                        for line in lp:
                            if line.startswith("ATOM"):
                                x = "%4.3f" % (float(line[30:37]) + shxb)
                                a = "%8s" % x

                                y = "%4.3f" % (float(line[38:45]) + shyb)
                                b = "%8s" % y
                                z = "%4.3f" % (float(line[46:53]) + shzb)
                                c = "%8s" % z
                                print b
                                line = line = line[0:29] + a + line[38:]
                                line = line = line[0:38] + b + line[46:]
                                line = line = line[0:46] + c + line[54:]

                        
                                of.write( line )
                        of.write("END")


class SpiderAnalyse(PyTool):
    args = [
    _( "dataset_dir", type="dir" ),
    _( "homologue_list", type="file" ),
    ]
    out = [
     _( "stat_file", file="statistics.csv" ),
     _( "res_file", file="overall.txt" )
    ]
    def func (self):
        fi=sorted((os.listdir(self.dataset_dir)), key= lambda x: int(x.split('_')[-1]))
        folders=os.listdir(self.dataset_dir)
        z=2
        homo = [line.rstrip() for line in open(self.homologue_list)]

        with open ("analysis.csv",'w') as ana:
            bestrmsds=[]
            bestrmsdsnh=[] #bestrmsdsno homologues
            wodens=[]
            wodensnh=[] #without density no homologues
            
            optimalrmsd=[]
            hc=0
            hc2=0
            print len(fi)
            for x in fi:
                y= int(filter(lambda x: x.isdigit(),x[-4:]))
                if y<999:
                    name="%s.%s" % (x,'pdb')
                    ori=os.path.join(self.dataset_dir,x,name)
                    orinpdb=numpdb.NumPdb(ori, {
                                    "phi_psi": False,
                                    "sstruc": False,
                                    "backbone_only": True,
                                    "protein_only": False,
                                    "detect_incomplete": False,
                                    "configuration": False,
                                    "info": False
                                    })
                    result={}
                    resultnh={}
                    #resultnh={}
                    resultlist=[]
                    
                    loops=os.path.join(self.dataset_dir,x,'loop_correl/crosscorrelation/loops')
                    ccsort=os.path.join(self.dataset_dir,x,'loop_correl/crosscorrelation/ccsort.cpv')
                    linkertxt=os.path.join(self.dataset_dir,x,'link_it/edited_linker.txt')
                    if os.path.isfile(ccsort):
                        with open (ccsort,'r') as hurz:
                            lines=hurz.readlines()
                            nr=lines[1].split()[0]
                            
                            #hurz=','.join(lines[1].split())
                            #print nr
                           # break
                        with open (linkertxt, 'r') as linki:
                            lineslinker=linki.readlines()
                        
                        h=os.listdir(loops)
                        for i in h:
                            if i.endswith('.pdb'):
                                #print i 
                                lpath=os.path.join(loops,i)
                                npdb = numpdb.NumPdb(lpath,{
                                            "phi_psi": False,
                                            "sstruc": False,
                                            "backbone_only": True,
                                            "protein_only": False,
                                            "detect_incomplete": False,
                                            "configuration": False,
                                            "info": False
                                            })
                                #print orinpdb['xyz']
                               # print lpath
                                lnr=int(i[0:3]) *4+1
                                #print orinpdb['xyz']
                                #print npdb['xyz']
                                rootm2=rmsd(orinpdb['xyz'],npdb['xyz'])
                                #print rootm2
                                result[i]=[rootm2]
                                result[i].append(lineslinker[lnr].strip())
                                resultlist.append(rootm2)
                                
                            else:
                                continue
                            
                        #den result dictonary auf eintraege resuzieren die nicht auf der Liste stehen:
                        resultnh=result.copy()
                        for key in result:
                            if result[key][1].upper() in homo:
                                del resultnh[key]
                        #print len(resultnh)
                        pdbkey="%s_%s.%s" % (nr.zfill(3),'bb','pdb')

                        
                        #abchecken ob 
                        #if result[pdbkey][1].upper() in homo:
                        #    hc+=1
                        #    print 'yes'
                        #else:
                        #    continue

                        best=sorted(resultlist)
                        #print result
                        print "\n",x,"\n"
                        #rmsd vonn loop nr 1
                        print "without desnity:" ,result['001_bb.pdb']
                        #rmsd vom hoechst geranketem loop
                        print "best:",pdbkey, result[pdbkey]
                        #der optimale rmsd: also einfach der kleinste aller 100
                        print "optimal",min(result, key=result.get),best[0]
                        nr_opt=min(result, key=result.get)
                        
                        i=1
                        
                        # an welcher Stelle is loop nr 1 gerankt
                        
                        for c,line in enumerate (lines, 1):
                            gu=line.split()
                            if gu[0]==str(nr_opt)[1:3]:
                                z=c
                                print "was ranked on number:", c-1
                            i+=1
                        #der schlechstest moegliche rmsd, also der letzte in der sortierten liste
                        print "worst",best[-1]
                        if pdbkey==min(result, key=result.get):
                            opti=1
                        else:
                            opti=0
                        if result[pdbkey][0]<result['001_bb.pdb'][0]:
                            better=1
                        else:
                            better=0
                        rline="%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % (x,result['001_bb.pdb'][0],pdbkey,result[pdbkey][0],min(result, key=result.get),best[0],z-1,best[-1],opti,better,'\n')
                        bestrmsds.append(result[pdbkey][0])
                        wodens.append(result['001_bb.pdb'][0])
                        optimalrmsd.append(best[0])
                        ana.write(rline)
                        
                        
                        ##############
                        
                        #schleife ueber die loops
                        
                        for y in range (1,100,1):
                            loop="%s_%s.%s" % (str(y).zfill(3),'bb','pdb')
                            print loop
                            if loop in resultnh:
                                wodensnh.append(result[loop][0])
                                break
                            else:
                                continue
                            
                        # das ganze fuer mit density
                        
                        for line in lines[1:]:
                            #lines.next()
                            no=line.split()[0]
                            lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                            
                            print 'richtig so?', result[lonam][1].upper()
                            if result[lonam][1].upper() not in homo:
                                print 'bloo'
                                bestrmsdsnh.append(result[lonam][0])
                                break
                            else:
                                continue
                            
        print hc
        print 'mean',sum(bestrmsds)/len(bestrmsds),'anzahl:',len(bestrmsds)
        print 'mean wo dens', sum(wodens)/len(wodens),'anzahl:',len(wodens)
        print 'best possible', sum (optimalrmsd)/len(optimalrmsd)
        print 'mean no homologues',sum(bestrmsdsnh)/len(bestrmsdsnh),'anzahl:',len(bestrmsdsnh)
        print 'mean wo dens no homologue', sum(wodensnh)/len(wodensnh),'anzahl:',len(wodensnh)

class LoopCrosscorrel2( PyTool ):
    args = [
        _( "mrc_file", type="file", ext="mrc" ),
        _( "pdb_file", type="file", ext="pdb" ),
        _( "loop_file", type="file", ext="pdb" ),
        _( "linkerinfo", type="file", ext="txt" ),
        _( "res1", type="sele" ),
        _( "res2", type="sele" ),
        _( "length", type="int", range=[1, 30] ),
        _( "resolution", type="float", range=[1, 10], step=0.1 ),
        _( "max_loops", type="int", range=[0, 500], default=100 )
    ]
    out = [
        _( "cropped_pdb", file="cropped.pdb" ),
        _( "ori_pdb_linker_file3", file="ori_pdb_linker_file3.pdb")
    ]
    tmpl_dir = TMPL_DIR

    def _init( self, *args, **kwargs ):
        self.spider_shift = SpiderShift(
            self.mrc_file, self.pdb_file,
            **copy_dict(
                kwargs, run=False, output_dir=self.subdir("shift")
            )
        )
        self.spider_convert = SpiderConvert(
            self.spider_shift.map_shift,
            **copy_dict(
                kwargs, run=False, output_dir=self.subdir("convert")
            )
        )
        self.spider_box = SpiderBox(
            self.spider_shift.map_shift,
            self.spider_convert.map_file,
            self.spider_shift.edited_pdb_file, self.res1, self.res2,
            self.length, self.resolution,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("box") )
        )
        self.pdb_box = SpiderPdbBox(
            self.spider_shift.edited_pdb_file,
            self.spider_box.box_file,
            self.spider_shift.map_shift,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("pdbbox"))
        )
        self.spider_delete_filled_densities = SpiderDeleteFilledDensities(
            self.spider_shift.map_shift,
            self.spider_box.box_map_file,
            self.pdb_box.edited_pdb_file,
            self.spider_box.box_file,
            self.resolution, self.res1, self.res2,
            **copy_dict(
                kwargs, run=False,
                output_dir=self.subdir("delete_filled_densities")
            )
        )
        self.spider_reconvert = SpiderReConvert(
            self.spider_box.box_file,
            self.spider_convert.map_file,
            self.spider_box.box_map_file,
            self.spider_convert.map_file,
            **copy_dict(
                kwargs, run=False, output_dir=self.subdir("reconvert")
            )
        )
        self.spider_crosscorrelation = SpiderCrosscorrelation(
            self.spider_convert.map_file,
            self.spider_delete_filled_densities.empty_map_file,
            self.spider_box.box_file,
            self.loop_file,
            self.linkerinfo,
            **copy_dict(
                kwargs, run=False, output_dir=self.subdir("crosscorrelation"),
                max_loops=self.max_loops
            )
        )
        self.output_files.extend( list( itertools.chain(
            self.spider_convert.output_files,
            self.spider_delete_filled_densities.output_files,
            self.spider_box.output_files,
            self.spider_reconvert.output_files,
            self.spider_crosscorrelation.output_files
        )))

    def func( self ):
        self._crop_pdb()
        self.spider_shift()
        self.spider_convert()
        self.spider_box()
        self.pdb_box()
        self.spider_delete_filled_densities()
        self.spider_reconvert()
        self.spider_crosscorrelation()

    def _post_exec( self ):
        self.backshift_linker()

    def _crop_pdb( self ):
        npdb = NumPdb( self.pdb_file )
        npdb.write(
            self.cropped_pdb,
            chain=self.res1["chain"],
            resno=[ self.res1["resno"] + 1, self.res2["resno"] - 1 ],
            invert=True
        )

    def backshift_linker( self ):
        shiftpath = self.relpath(self.spider_shift.output_dir)
        shiftfile = "%s/%s" % (shiftpath, 'shift.cpv')
        with open( shiftfile, "r" ) as rs:
            file_lines = rs.readlines()
            loc_file = file_lines[1].strip()
            shxb = float(loc_file[5:13]) * -1
            shyb = float(loc_file[19:27]) * -1
            shzb = float(loc_file[33:42]) * -1

        loop_dir = self.relpath(self.spider_crosscorrelation.loop_dir)
        outputdir = self.subdir("oriloops")
        with open( self.linkerinfo, 'r') as li:
            lilines = li.readlines()
        with open( self.loop_file, 'r' ) as lf:
            with open( self.ori_pdb_linker_file3, 'w' ) as fp_out:
                for line in lf:
                    if line.startswith("ATOM"):
                        x = "%4.3f" % (float(line[30:37]) + shxb)
                        a = "%8s" % x

                        y = "%4.3f" % (float(line[38:45]) + shyb)
                        b = "%8s" % y
                        z = "%4.3f" % (float(line[46:53]) + shzb)
                        c = "%8s" % z
                        print b
                        line = line = line[0:29] + a + line[38:]
                        line = line = line[0:38] + b + line[46:]
                        line = line = line[0:46] + c + line[54:]


                        fp_out.write( line )
                    else:
                        continue
                fp_out.write( "END" )
        for i in os.listdir(loop_dir):
            if i.endswith(".pdb"):
                loopfile = "%s/%s" % (loop_dir, i)

                outfile = "%s/%s" % (outputdir, i)
                with open( loopfile, 'r' ) as lp:
                    lonr = int(i[0:3])
                    lino = (lonr * 4) - 2

                    title = "%s     %s %s ,%s ,%s ,%s ,%s" % (
                        "TITLE", "linker", lonr, lilines[lino].strip(),
                        lilines[lino + 1].strip(), lilines[lino + 2].strip(),
                        lilines[lino + 3])

                    with open( outfile, 'w' ) as of:
                        of.write(title)
                        for line in lp:
                            if line.startswith("ATOM"):
                                x = "%4.3f" % (float(line[30:37]) + shxb)
                                a = "%8s" % x

                                y = "%4.3f" % (float(line[38:45]) + shyb)
                                b = "%8s" % y
                                z = "%4.3f" % (float(line[46:53]) + shzb)
                                c = "%8s" % z
                                print b
                                line = line = line[0:29] + a + line[38:]
                                line = line = line[0:38] + b + line[46:]
                                line = line = line[0:46] + c + line[54:]


                                of.write( line )
                        of.write("END")


class LoopCrosscorrel3( PyTool ):
    args = [
        _( "mrc_file", type="file", ext="mrc" ),
        _( "pdb_file", type="file", ext="pdb" ),
        _( "loop_file", type="file", ext="pdb" ),
        _( "linkerinfo", type="file", ext="txt" ),
        _( "res1", type="sele" ),
        _( "res2", type="sele" ),
        _( "length", type="int", range=[1, 30] ),
        _( "resolution", type="float", range=[1, 10], step=0.1 ),
        _( "max_loops", type="int", range=[0, 500], default=100 )
    ]
    out = [
        _( "cropped_pdb", file="cropped.pdb" ),
        _( "ori_pdb_linker_file3", file="ori_pdb_linker_file3.pdb")
    ]
    tmpl_dir = TMPL_DIR
    def _init( self, *args, **kwargs ):
        #self.spider_shift = SpiderShift(
        #    self.mrc_file,self.pdb_file,
        #    **copy_dict( 
        #        kwargs, run=False, output_dir=self.subdir("shift") 
        #    )
        #)
        #self.spider_convert = SpiderConvert( 
        #    self.spider_shift.map_shift, 
        #    **copy_dict( 
        #        kwargs, run=False, output_dir=self.subdir("convert") 
        #    )
       # )
        self.spider_box = SpiderBox(
            self.mrc_file,
            #self.spider_shift.map_shift,
            "/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/20/files/mapupload.cpv",
            #self.spider_convert.map_file, 
            #self.spider_shift.edited_pdb_file
            "/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/20/files/edited.pdb", self.res1, self.res2, 
            self.length, self.resolution,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("box") )
        )
        self.pdb_box =  SpiderPdbBox(
           #self.spider_shift.edited_pdb_file,
           "/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/20/files/edited.pdb",
           self.spider_box.box_file,
           self.mrc_file,
           #self.spider_shift.map_shift,
           **copy_dict( kwargs, run=False, output_dir=self.subdir("pdbbox"))
        )
        self.spider_delete_filled_densities = SpiderDeleteFilledDensities(
            self.mrc_file,self.spider_box.box_map_file,self.pdb_box.edited_pdb_file,self.spider_box.box_file,
            
            #self.spider_shift.map_shift,self.spider_box.box_map_file, self.pdb_box.edited_pdb_file, self.spider_box.box_file,
            self.resolution,self.res1,self.res2,
            **copy_dict( 
                kwargs, run=False, 
                output_dir=self.subdir("delete_filled_densities") 
            )
        )
        self.spider_reconvert = SpiderReConvert(
            self.spider_box.box_file,
            "/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/20/files/mapupload.cpv",
            #self.spider_convert.map_file,
            self.spider_box.box_map_file,
            "/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/20/files/mapupload.cpv",
            #self.spider_convert.map_file,
            **copy_dict( 
                kwargs, run=False, output_dir=self.subdir("reconvert") 
            )
        )
        self.spider_crosscorrelation = SpiderCrosscorrelation(
            "/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/20/files/mapupload.cpv",
            #self.spider_convert.map_file, 
            self.spider_delete_filled_densities.empty_map_file, 
            self.spider_box.box_file, 
            self.loop_file,
            self.linkerinfo,
            **copy_dict( 
                kwargs, run=False, output_dir=self.subdir("crosscorrelation"),
                max_loops=self.max_loops
            )
        )
        self.output_files.extend( list( itertools.chain(
            #self.spider_convert.output_files, 
            self.spider_delete_filled_densities.output_files,
            self.spider_box.output_files,
            self.spider_reconvert.output_files,
            self.spider_crosscorrelation.output_files
        )))
    def func( self ):
        self._crop_pdb()
        #self.spider_shift()
        #self.spider_convert()
        self.spider_box()
        self.pdb_box ()
        self.spider_delete_filled_densities()
        self.spider_reconvert()
        self.spider_crosscorrelation()
    
    def _post_exec( self ):
        self.backshift_linker()
        
    def _crop_pdb( self ):
        npdb = NumPdb( self.pdb_file )
        npdb.write( 
            self.cropped_pdb, 
            chain=self.res1["chain"], 
            resno=[ self.res1["resno"]+1, self.res2["resno"]-1 ],
            invert=True
        )
    
    def backshift_linker ( self ) :
        #shiftpath=self.relpath(self.spider_shift.output_dir)
        shiftfile='/home/jochen/work/fragfit/benchmark/sim/gpcrcomplex/20/files/shift.cpv'
        #shiftfile="%s/%s" % (shiftpath,'shift.cpv')
        with open (shiftfile, "r") as rs:
            file_lines=rs.readlines()
            loc_file = file_lines[1].strip()
            shxb=float(loc_file[5:13])*-1
            shyb=float(loc_file[19:27])*-1
            shzb=float(loc_file[33:42])*-1

            backshift= np.array([shxb,shyb,shzb])
   
        loop_dir=self.relpath(self.spider_crosscorrelation.loop_dir)
        outputdir=self.subdir("oriloops")
        with open ( self.linkerinfo , 'r') as li:
            lilines=li.readlines()
        with open (self.loop_file, 'r') as lf:
            with open (self.ori_pdb_linker_file3, 'w') as fp_out:
                for line in lf:
                    if line.startswith("ATOM"):
                        x = "%4.3f" % (float(line[30:37]) + shxb)
                        a = "%8s" % x

                        y = "%4.3f" % (float(line[38:45]) + shyb)
                        b = "%8s" % y
                        z = "%4.3f" % (float(line[46:53]) + shzb)
                        c = "%8s" % z
                        print b
                        line = line = line[0:29] + a + line[38:]
                        line = line = line[0:38] + b + line[46:]
                        line = line = line[0:46] + c + line[54:]

        
                        fp_out.write( line )
                    else:
                        continue
                fp_out.write( "END" )
        for i in os.listdir(loop_dir):
            if i.endswith(".pdb"):
                loopfile="%s/%s" % (loop_dir,i)
                
                outfile="%s/%s" %(outputdir,i)
                with open (loopfile, 'r') as lp:
                    lonr=int(i[0:3])
                    lino=(lonr*4)-2

                    title= "%s     %s %s ,%s ,%s ,%s ,%s" % ("TITLE","linker",lonr,lilines[lino].strip(),lilines[lino+1].strip(),lilines[lino+2].strip(),lilines[lino+3])

                    with open (outfile, 'w') as of:
                        of.write(title)
                        for line in lp:
                            if line.startswith("ATOM"):
                                x = "%4.3f" % (float(line[30:37]) + shxb)
                                a = "%8s" % x

                                y = "%4.3f" % (float(line[38:45]) + shyb)
                                b = "%8s" % y
                                z = "%4.3f" % (float(line[46:53]) + shzb)
                                c = "%8s" % z
                                print b
                                line = line = line[0:29] + a + line[38:]
                                line = line = line[0:38] + b + line[46:]
                                line = line = line[0:46] + c + line[54:]

                                of.write( line )
                        of.write("END")
class SpiderAnalyse2(PyTool):
    args = [
    _( "dataset_dir", type="dir" ),
    _( "loop_length", type="int" ),
    _( "homologue_list", type="file" ),
    _( "header_pdb", type="file" ),
    _("analysis_flag", type="int")
    ]
    out = [
     _( "stat_file", file="statistics.csv" ),
     _( "res_file", file="overall.txt" )
    ]
    def func (self):
        helix={}
        sheet={}
        fullsheet={}
        sheetnames=[]
        with open (self.header_pdb) as header:
            headerlines=header.readlines()
            xyz=0
            for lines in headerlines:
                if lines.startswith  ("HELIX"):
                    ident="%s_%s_%s" % ('HELIX',lines[19:20],lines[8:10])
                    helix[ident]=[int(lines[72:76])]
                    helix[ident].append(lines[8:10])
                    helix[ident].append(lines[31:32])
                    helix[ident].append(int(lines[22:25]))
                    helix[ident].append(lines[34:37])
                if lines.startswith  ("SHEET"):
                    #print int(lines[34:37])
                    xyz+=1
                    #hurzi=lines.split()
                    #print hurzi[4]
                    ident=lines[12:14],lines[8:10]#"%s_%s_%s.%s" % ("SHEET",lines[21:22],lines[8:10],xyz)
                    sheetlen=(int(lines[34:37])-int(lines[23:26]))+1
                    #print sheetlen
                    sheet[ident]=[sheetlen]#[int(lines[72:76])]
                    sheet[ident].append(lines[8:10])
                    sheet[ident].append(lines[21:22])
                    sheet[ident].append(int(lines[23:26]))
                    sheet[ident].append(int(lines[34:37]))
                    sheet[ident].append(lines[12:14])
                    sheetnames.append(lines[12:14])
        for hurz in sorted(set(sheetnames)):
            go={}
            for key in sorted(sheet.iterkeys()):
                
                if key[0]==hurz:

                    go[key]=sheet[key]

            for i in range (0,len(go)-1):

                if go.values()[i] [2]==go.values()[i+1] [2] and abs(go.values()[i] [4]-go.values()[i+1] [4])<10 :
 
                    sheetlist=[]
                    sheetlist.extend([go.values()[i] [4],go.values()[i] [3],go.values()[i+1] [3],go.values()[i+1] [4]])
                    identfull= "%s_%s_%s" % (hurz,go.values()[i] [2],min(sheetlist))
                    fullsheet[identfull]=[go.values()[i] [2]]
                    fullsheet[identfull].append(min(sheetlist))
                    fullsheet[identfull].append(max(sheetlist))
                    fullsheet[identfull].append(max(sheetlist)-min(sheetlist))
        for key in sorted(fullsheet):
          #  print key, fullsheet[key]
          continue
        with open ('loopana_B.csv','w') as analyse, open ('values2.csv', 'w') as valu, open ('cases.txt', 'w') as cases, open ('bettercc.txt', 'w') as bettercc, open ('diff.csv','w') as diffile:
            title="%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % ("looplength","fragfit","std rmsd with density",
                                                                   "rmsd linkit","std rmsd without density",
                                                                   "rmsd fragfit no homologues","std fragfit no homologues",
                                                                   "rmsd linkit no homologs","std rmsd linkit no homologs",
                                                                   "rmsd fragfit no homologues 75","rmsd fragfit no homologues 90","std fragfit no linker homologues",
                                                                   "rmsd linkit no linker homologs","std linkit no linker homologs",
                                                                   "Anzahl Suchen","optimal","mean_correl","extragut","beste five","best ten","linkit_helical",
                                                                   "linkit_other","mean rank","median fragfit","tms scrore","Fragfit 50","Fragfit 20",
                                                                   'ori cc mean','Diff to oricc','diff to oricc homo75','diff 75','ccbestfive','meanccbest','fragfitheli','fragfitother','fragfitbeta','sequenceidentity no homo','sequence identity full','fragfithel ori','fragfitbeta ori','fragfitother ori','anteil perfekt','anteilperfekt nh','linkitbetter','fragfitbetter','linkitfirst','fragfitfirst',
                                                                   'linkit>5','fragfitof>5','linkitbad','fragfitbad','conuter bad',
                                                                   "helix_link_it","helix_fragfit",'anzahl helix',"sheet_linkit",'Sheet_fragfit','anzahl sheet','fullsheetlinkit','fullsheetfragfit','anzahl fullsheet','tmdo.5','tmdo_homo.5','seqidff','seqIDFsNh','seqIDoptimal','cc75','\n')
            analyse.write(title)
            rmsdval={
            '5': 1.7857142857,
'6':2.1428571429,
'7':2.5,
'8':2.8571428571,
'9':3.2142857143,
'10':3.5714285714,
'11':3.9285714286,
'12':4.2857142857,
'13':4.6428571429,
'14':5,
'15':5.3571428571,
'16':5.7142857143,
'17':6.0714285714,
'18':6.4285714286,
'19':6.7857142857,
'20':7.1428571429,
'21':7.5,
'22':7.8571428571,
'23':8.2142857143,
'24':8.5714285714,
'25':8.9285714286,
'26':9.2857142857,
'27':9.6428571429,
'28':10,
'29':10.3571428571,
'30':10.7142857143,
'31':11.0714285714,
'32':11.4285714286,
'33':11.7857142857,
'34':12.1428571429
  
            }
            correlation={
            '5':0.568783917098,
            '6':0.56042945098,
            '7':0.523105761658,
            '8':0.463180293532,
            '9':0.41430411,
            '10':0.367170015707,
            '11':0.329382698022,
            '12':0.298228513661,
            '13':0.250647477612,
            '14':0.222953401736,
            '15':0.199863570246,
            '16':0.174128295472,
            '17':0.167516662037,
            '18':0.10233948813,
            '19':0.0990273932505,
            '20':0.0832893506279,
            '21':0.0721325818392,
            '22':0.0700303510774,
            '23':0.0653087346573,
            '24':0.0603087346573,
            '25':0.05819717193,
            '26':0.0612093805726,
            '27':0.0621014396482,
            '28':0.0568244603245,
            '29':0.0757444702719,
            '30':0.0770906372842,
            '31':0.0669438991446,
            '32':0.0775876085589,
            '33':0.0793356798089,
            '34':0.0514788622987,
            '35':0.0514788622987
            
            }
            for l in range(5,35,1):
                indifile="%s_%s.%s" % (l,'file','csv')
                t_test="%s_%s.%s" % (l,'t_test','txt')
                with open (indifile, 'w') as loofile:
                    title_loo="%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % ('chain','name','taken','ori_pdb','seq_identit','rmsd','SSE flag','diffccabs','tmscore','diffreal','bestes_modell','bestrmsd','fragfitbest','fragfitbesrmsd','bestpossibleall?','bestlinkitnh','bestpossiblenh','linkitbest','fragfitbest','linkit>5','fragfitof>5','linkitbad','fragfitbad','\n')
                    loofile.write(title_loo)
                    diffile.write(str(l) + ',')
                    scorr=correlation[str(l)] #standartdcorreltion fuer die looplaenge
                    corrband1=scorr+scorr/5
                    corrband2=scorr-scorr/5
                    #print corrband2, corrband1
                    #break
                    #corrrange=range(corrband2,corrband1)
                    chainfolders=os.listdir(self.dataset_dir)
                    bestrmsds=[]
                    bestrmsdsnh=[] #bestrmsdsno homologues
                    wodens=[]
                    wodensnh=[] #without density no homologues
                    wodensnhli=[]
                    bestrmsdnhli=[] #no homologue linker
                    tmsnhli=[]
                    meancorrel=[]
                    ori_meancorrel=[]
                    optimalrmsd=[]
                    extragut=[]
                    bestfiveall=[]
                    besttenall=[]
                    linkitheli=[]
                    linkitother=[]
                    fragfithomo90=[]
                    rank=[]
                    result50=[]
                    result20=[]
                    hc=0
                    hc2=0
                    sucounter=0
                    sucounternh=0
                    greatdiff=[]
                    or_diffcc=[]
                    ordiff75=[]
                    ordiffh75full=[]
                    diffbestfive=[]
                    meancc75=[]
                    fragfitheli=[]
                    fragfitother=[]
                    fragfitbeta=[]
                    sequenceidentity=[]
                    seqidentfull=[]
                    seqidfffull=[]
                    fragfithelior=[]
                    fragfitbetaor=[]
                    fragitotheror=[]
                    linkitonlybetter=[]
                    fragfitonlybetter=[]
                    fragfitfirst=[]
                    linkitfirst=[]
                    fragfitfirstbad=[]
                    linkitfirstbad=[]
                    fragfitfirstbadf=[]
                    linkitfirstbadf=[]
                    helixlinkit=[]
                    helixfragfit=[]
                    sheetlinkit=[]
                    sheetfragfit=[]
                    fullsheetlinkit=[]
                    fullsheetfragfit=[]
                    tmscounter=[]
                    tmscounterho=[]
                    fsnohomoid=[]
                    seqidoptimal=[]
                    cc75=[]
                    badf=0

                    for g in chainfolders:
                        #g='A'
                        #print g
                        loopfolder=os.path.join(self.dataset_dir,g,str(l))
                        if os.path.isdir(loopfolder):
                            fi=sorted((os.listdir(loopfolder)), key= lambda x: int(x.split('_')[-1]))
                            #folders=os.listdir(self.dataset_dir)
                            z=2
                            homo = [line.rstrip() for line in open(self.homologue_list)]
                           # print loopfolder
                            with open ("analysis.csv",'w') as ana:
                                
                                #print fi
                                for x in fi:
                                   # print x
                                    y= int(filter(lambda x: x.isdigit(),x[-4:]))
                                    y1=int(filter(lambda x: x.isdigit(),x[:3]))
                                   #print y1,y
                                    if y<999:
    
                                        name="%s.%s" % (x,'pdb')
                                        namezipped="%s.%s" % (x,'pdb.zip')
                                        zout=os.path.join(loopfolder,x)
                                        ori=os.path.join(loopfolder,x,name)
                                        orizipped=os.path.join(loopfolder,x,namezipped)
                                        if os.path.isfile(orizipped):
                                            #print "hallo"
                                            with zipfile.ZipFile(orizipped,'r') as zo:
                                                zo.extractall(loopfolder)
                                            #PdbUnzip(orizipped)
                                        folo=os.path.join(loopfolder,x)
                                        #print ori
                                        
                                        orinpdb=numpdb.NumPdb(ori, {
                                                        "phi_psi": False,
                                                        "sstruc": False,
                                                        "backbone_only": True,
                                                        "protein_only": False,
                                                        "detect_incomplete": False,
                                                        "configuration": False#,
                                                        #"info": False
                                                        })
                                        ccdic={}
                                        result={}
                                        resultnh={}
                                        #resultnh={}
                                        resultlist=[]
                                        seq1=orinpdb.sequence()
                                        loops=os.path.join(loopfolder,x,'loop_correl/crosscorrelation/loops_repair')
                                        ccsort=os.path.join(loopfolder,x,'loop_correl/crosscorrelation/ccsort.cpv')
                                        linkertxt=os.path.join(loopfolder,x,'link_it/edited_linker.txt')
                                        oricc=os.path.join(loopfolder,x,'oricrosscorrelation.cpv')
                                        if os.path.isfile(ccsort):
                                            Dssp(ori,output_dir=folo)
                                            oridssp='%s%s' % (ori[:-3],'dssp')
                                            orirec=parse_dssp(oridssp)
                                            with open (ccsort,'r') as hurz:
                                                lines=hurz.readlines()
                                               # print ori
                                                nr=lines[1].split()[0]
                                                corr=lines[1].split()[2]
                                                meancorrel.append(float(corr))
                                                for j in lines[1:]:
                                                    u=int(j.split()[0])
                                                    
                                                    cci=j.split()[2]
                                                    
                                                    ccdic[u]=float(cci)
                                               
                                            with open (linkertxt, 'r') as linki:
                                                lineslinker=linki.readlines()
                                            #with open (oricc,'r') as oriloopcc:
                                            #    linesoricc=oriloopcc.readlines()
                                            #    oriccval=linesoricc[1].split()[2]
                                            #    diffcc=abs(float(corr)-float(oriccval))
                                            #    if float(corr) > float(oriccval):
                                            #        dro=  "%s/%s %s %s%s" % (loopfolder,x,float(corr), float(oriccval), '\n')
                                            #        bettercc.write( dro )
                                            #    ori_meancorrel.append (float(oriccval))
                                            #    or_diffcc.append(diffcc)
                                            #    #print oriccval
                                            oriccval=0
                                            h=os.listdir(loops)
                                            for i in h:
                                                if i.endswith ('.zip'):
                                                    pdbzipname=os.path.join(loops,i)
                                                    #print pdbzipname
                                                    with zipfile.ZipFile(pdbzipname,'r') as zi:
                                                        zi.extractall(loopfolder)
                                                    os.remove(pdbzipname)
                                            h=os.listdir(loops)
                                            for i in h:
                                                if i.endswith('.pdb'):
                                                    #print i
                                                    lino=int(i[0:3])
                                                    lpath=os.path.join(loops,i)
                                                    npdb = numpdb.NumPdb(lpath,{
                                                                "phi_psi": False,
                                                                "sstruc": False,
                                                                "backbone_only": True,
                                                                "protein_only": False,
                                                                "detect_incomplete": False,
                                                                "configuration": False,
                                                                "info": False
                                                                })
                                                    #print npdb.get
                                                    res1=int(x.split('_')[0])
                                                    res2=int(x.split('_')[1])
                                                   # print res1, res2
                                                    npdbshort=npdb.copy(resno=[res1,res2])
                                                    lnr=int(i[0:3]) *4+1
                                                    seqnr=int(i[0:3]) *4
    
                                                    rootm2=rmsd(orinpdb['xyz'],npdbshort['xyz'])
    
                                                    try:
    
                                                        tms=tmscore(orinpdb['xyz'],npdbshort['xyz'])
                                                    except:
                                                        tms=0
                                                        #print l,x,i
                                                      
                                                    seq2=lineslinker[seqnr].strip()
                                                    
                                                    #print seq1
                                                    #sequence identity
                                                    si=0
                                                    for sf in range(0, len(seq1),1):
                                                        if seq1[sf]==seq2[sf+1]:
                                                            si+=1
                                                    sqi=(si/len(seq1))*100
                                                    #print seq1,"and",seq2[1:-1]
                                                    print sqi
                                                        
                                                    
                                                        
                                                    #print rootm2
                                                    result[i]=[rootm2]
                                                    result[i].append(lineslinker[lnr].strip())
                                                    result[i].append(sqi)
                                                    result[i].append(tms)
                                                    result[i].append(ccdic[lino])
                                                    #result[i].append(oriccval)
                                                    resultlist.append(rootm2)
                                                    
                                                else:
                                                    continue
                                                
                                            #den result dictonary auf eintraege resuzieren die nicht auf der Liste stehen:
                                            resultnh=result.copy()
                                            for key in result:
                                                if result[key][1].upper() in homo:
                                                    del resultnh[key]
                                            #einen dictonary erstellen wo die homologen linker nicht drin sind:
                                            resultnhli=result.copy()
                                            for key in result:
                                                if result[key][2] >=75:
                                                    del resultnhli[key]
                                    
                                            pdbkey="%s_%s.%s" % (nr.zfill(3),'bb','pdb')
    
                    
                                            best=sorted(resultlist)
    
                                            nr_opt=min(result, key=result.get)
                                            
                                            i=1
                                            
                                            # an welcher Stelle is loop nr 1 gerankt
                                            
                                            for c,line in enumerate (lines, 1):
                                                gu=line.split()
                                                if gu[0]==str(nr_opt)[1:3]:
                                                    z=c
                                                    rank.append(c-1)
    
                                                i+=1
    
                                            if pdbkey==min(result, key=result.get):
                                                opti=1
                                            else:
                                                opti=0
                                            if result[pdbkey][0]<result['001_bb.pdb'][0]:
                                                better=1
                                            else:
                                                better=0
                                            if result[pdbkey][0]-result['001_bb.pdb'][0]>4:
                                                cases.write(ori + os.linesep)
                                                #print ori
                                            #print 'min', min(result, key=result.get)
                                            #print 'result',result['001_bb.pdb'][0]
                                            #print 'pdbkey' ,pdbkey
                                            rline="%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % (x,result['001_bb.pdb'][0],pdbkey,result[pdbkey][0],min(result, key=result.get),best[0],z-1,best[-1],opti,better,'\n')
                                            bestrmsds.append(result[pdbkey][0])
                                            seqidentfull.append(result[pdbkey][2])
                                            seqidfffull.append(result['001_bb.pdb'][2])
                                            wodens.append(result['001_bb.pdb'][0])
                                            optimalrmsd.append(best[0])
                                            seqidoptimal.append(result[nr_opt][2])
                                            
                                            ###nur extragute
                                            if corrband2<float(corr) <corrband1 and result[pdbkey][2]<75:
                                                extragut.append(result[pdbkey][0])
                                            #
                                            ana.write(rline)
                                            
                                            
                                            ##############
                                            
                                            #schleife ueber die loops
                                            
                                            for y in range (1,100,1):
                                                loop="%s_%s.%s" % (str(y).zfill(3),'bb','pdb')
                                               # print loop
                                                if loop in resultnh:
                                                    wodensnh.append(result[loop][0])
                                                    break
                                                else:
                                                    continue
                                                
                                            # das ganze fuer mit density
                                            
                                            for line in lines[1:]:
                                                #lines.next()
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                
                                               # print 'richtig so?', result[lonam][1].upper()
                                                if result[lonam][1].upper() not in homo:
                                                    #print 'bloo'
                                                    bestrmsdsnh.append(result[lonam][0])
                                                    break
                                                else:
                                                    continue
                                            #homologe linker aussortieren:
                                            
                                            for y in range (1,100,1):
                                                loop="%s_%s.%s" % (str(y).zfill(3),'bb','pdb')
                                                #print loop
                                               # print loop
                                                if loop in resultnhli:
                                                    bestlinkitnohomo=loop
                                                    wodensnhli.append(result[loop][0])
                                                    fsnohomoid.append(result[loop][2])
                                                    #print loop
                                                    hfj=os.path.join(loops,loop)
                                                    hjd='%s%s' % (hfj[:-3],'dssp')
                                                    Dssp(hfj,output_dir=loops)
                                                    dsrec=parse_dssp(hjd)
                                                    
                                                  #  print hjd
                                                    
                                                    #with open ('hjd') as ds:
                                                    #    gd=ds.readlines()
                                                    #    for line in gd:
                                                    #        if line.split() []
                                                    if any('H'==e[3]for e in dsrec[1:])==True:
                                                        #print 'found helical loop'
                                                        linkitheli.append(result[loop][0])
                                                        #fragfitheli.append
                                                    else:
                                                        linkitother.append(result[loop][0])
                                                    #for i in dsrec[1:]:
                                                    #    if i[3]=='H':
                                                    #        print i[3]
                                                    #print "dssp", dsrec[1]
                                                    break
                                                else:
                                                    continue
                                            for line in lines[1:]:
                                                #lines.next()
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                hfj2=os.path.join(loops,lonam)
                                                Dssp(hfj2,output_dir=loops)
                                                hjd2='%s%s' % (hfj2[:-3],'dssp')
                                                dsrec2=parse_dssp(hjd2)
                                               # print dsrec2
                                                ccsortbest=result[lonam][0]
                                                if result[lonam][2]<75:                                               
                                                    if any('H'==e[3]for e in dsrec2[1:])==True:
                                                        fragfitheli.append(result[lonam][0])
                                                    elif any ('E'==e[3]for e in dsrec2[1:])==True:
                                                        fragfitbeta.append(result[lonam][0])
                                                    else:
                                                        fragfitother.append(result[lonam][0])
                                                    break
                                                else:
                                                    continue
                                            optrmsdnhnr=min(resultnhli,key=resultnhli.get)
                                            optrmsdnhl=result[optrmsdnhnr][0]
                                           # print optrmsdnhl
                                            for line in lines[1:]:
                                                #lines.next()
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                if result[lonam][2]<75:
                                                    
                                                    bestrmsdnhli.append(result[lonam][0])
                                                    bestfragnoh=result[lonam][0]
                                                    tmsnhli.append(result[lonam][3])
                                                    cc75.append(result[lonam][4])
                                                    if result[lonam][3] >0.5:
                                                        tmscounterho.append(1)
                                                    else:
                                                        tmscounterho.append(0)
                                                    
                                                    meancc75.append(result[lonam][4])
                                                    diffh75= abs(float(ccdic[int(no)])-float(oriccval))
                                                    diffh75full=float(ccdic[int(no)])-float(oriccval)
                                                    ordiff75.append( diffh75 )
                                                    ordiffh75full.append(diffh75full )
                                                    suflagnh=0
                                                    if bestlinkitnohomo==optrmsdnhl:
                                                        continue
                                                    else:
                                                        linkitonlybetter.append(result[bestlinkitnohomo][0])
                                                        fragfitonlybetter.append(result[lonam][0])
                                                    if result[lonam][0]==optrmsdnhl:
                                                        sucounternh+=1
                                                        suflagnh=1
                                                    else:
                                                        suflagnh=0   
                                            
                                                
                                                    
                                                    try:
                                                        zhf=0
                                                        if any ('H'==e[3] for e in orirec[1:])==True:
                                                            fragfithelior.append(result[lonam][0])
                                                            zhf=1
                                                        elif any ('E'==e[3] for e in orirec[1:])==True:
                                                            fragfitbetaor.append(result[lonam][0])
                                                            zhr=2
                                                        else:
                                                            fragitotheror.append(result[lonam][0])
                                                            zhf=3
                                                        diffile.write(str(float(ccdic[int(no)])-float(oriccval)) +',')
                                                        sequenceidentity.append(result[lonam][2])
                                                        
                                                        #if math.isnan(diffh75):
                                                        #    print loopfolder
                                                        #    break
                                                        #print diffh75
                                                        break
                                                    except:
                                                        continue
                                                else:
                                                    continue
                                            suflag=0
                                            if result[pdbkey][0] ==best[0]:
                                                suflag=1
                                                sucounter +=1
                                            else:
                                                suflag=0
                                            linkitfirst.append(result['001_bb.pdb'][0])
                                            fragfitfirst.append(ccsortbest)
                                            if result['001_bb.pdb'][3] >0.5:
                                                tmscounter.append(1)
                                            else:
                                                tmscounter.append(0)
                                            if result['001_bb.pdb'][0]> 5:
                                                linkitfirstbad.append(result['001_bb.pdb'][0])
                                                fragfitfirstbad.append(ccsortbest)
                                            if result['001_bb.pdb'][0]>rmsdval[str(l)]:
                                                linkitfirstbadf.append(result['001_bb.pdb'][0])
                                                fragfitfirstbadf.append(ccsortbest)
                                                badf+=1
                                            looline= "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % (g,x,lonam,result[lonam][1],result[lonam][2],result[lonam][0],zhf,str(float(ccdic[int(no)])-float(oriccval)),result[lonam][3],diffh75full,nr_opt,best[0],pdbkey,result[pdbkey][0],suflag,optrmsdnhl,suflagnh,result['001_bb.pdb'][0],ccsortbest,'\n')
                                            loofile.write(looline)
                                            for line in lines[1:]:
                                                #lines.next()
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                
                                                if result[lonam][2]<90:
                                                    fragfithomo90.append(result[lonam][0])
                                                    #tmsnhli.append(result[lonam][3])
                                                    break
                                                else:
                                                    continue
                                            
                                            for line in lines[1:]:
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                if  int(no)<50 and result[lonam][2]<75:
                                                    result50.append(result[lonam][0])
                                                    break
                                                else:
                                                    continue
                                            
                                            for line in lines[1:]:
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                if  int(no)<20 and result[lonam][2]<75:
                                                    result20.append(result[lonam][0])
                                                    break
                                                else:
                                                    continue
                                            #die besten five alternative
                                            bestfive={}
                                            for line in lines[1:]:
                                            #    #lines.next()
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                if result[lonam][2]<75:
                                                    
                                                    bestfive[lonam]=(result[lonam])
                                                   
                                                    if len(bestfive)==5:
                                                        #print bestfive
                                                        
                                                       # print sorted(bestfive.values())[0] 
                                                        bestfiveall.append(sorted(bestfive.values())[0] [0])
                                                        diffbestfive.append(sorted(bestfive.values())[0] [4])
                                                        break
                                                else:
                                                    continue
                                            #die besten zehn
                                            #bestfive=[]
                                            #for line in lines[1:]:
                                            #    #lines.next()
                                            #    no=line.split()[0]
                                            #    lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                            #    if result[lonam][2]<75:
                                            #        bestfive.append(float(result[lonam][0]))
                                            #        if len(bestfive)==5:
                                            #            #print sorted(bestfive)
                                            #            bestfiveall.append(sorted(bestfive)[0])
                                            #            break
                                            #    else:
                                            #        continue
                                            #die besten zehn
                                            bestten=[]
                                            for line in lines[1:]:
                                                #lines.next()
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                if result[lonam][2]<75:
                                                    bestten.append(float(result[lonam][0]))
                                                    if len(bestten)==10:
                                                        #print sorted(bestfive)
                                                        besttenall.append(sorted(bestten)[0])
                                                        break
                                                else:
                                                    continue       # print '.'
                                            #helix_sheet_by header
                                            for key in helix:
                                                        if helix[key] [0] >=l and helix[key] [2]==g and helix[key] [3] in range (y1-2,y1+2):
                                                          #  print key, helix[key], g, y1
                                                            helixlinkit.append(result['001_bb.pdb'][0])
                                                            helixfragfit.append(ccsortbest)
                                                            
                                            for key in sheet:
                                                if sheet[key] [0] >=l and sheet[key] [2]==g and sheet[key] [3] in range (y1-2,y1+2):
                                                   # print key, sheet[key], g, y1
                                                    sheetlinkit.append(result['001_bb.pdb'][0])
                                                    sheetfragfit.append(ccsortbest)
                                            for key in fullsheet:
                                                if fullsheet[key] [3] >=l and fullsheet[key] [0]==g and fullsheet[key] [1] in range (y1-2,y1+2):
                                                    #print key, sheet[key], g, y1
                                                    fullsheetlinkit.append(result['001_bb.pdb'][0])
                                                    fullsheetfragfit.append(ccsortbest)
                if len(linkitheli)==0:
                    linkitheli.append(0)
                if len(extragut)==0:
                    extragut.append(0)
                if len(bestrmsds)==0:
                    bestrmsds.append(0)
                if  len(wodens)==0:
                    wodens.append(0)
                if len(optimalrmsd)==0:
                    optimalrmsd.append(0)
                if len(bestrmsdsnh)==0:
                    bestrmsdsnh.append(0)
                if len(wodensnh)==0:
                    wodensnh.append(0)
                if len(bestrmsdnhli)==0:
                    bestrmsdnhli.append(0)
                if len(wodensnhli)==0:
                    wodensnhli.append(0)
                if  len(meancorrel)==0:
                    meancorrel.append(0)
                if len(bestfiveall)==0:
                    bestfiveall.append(0)
                if len(linkitother)==0:
                    linkitother.append(0)
                if len(rank)==0:
                    rank.append(0)
                if len(tmsnhli)==0:
                    tmsnhli.append(0)
                if len(result50)==0:
                    result50.append(0)
                if len(fragfithomo90)==0:
                    fragfithomo90.append(0)
                if len (result20)==0:
                    result20.append(0)
                if len (besttenall)==0:
                    besttenall.append(0)
                if len (ori_meancorrel)==0:
                    ori_meancorrel.append(0)
                if len (or_diffcc)==0:
                    or_diffcc.append(0)
                if len (ordiff75)==0:
                  #  print loopfolder
                    ordiff75.append(0)
                if len ( ordiffh75full)==0:
                     ordiffh75full.append(0)
                if len (diffbestfive)==0:
                    diffbestfive.append(0)
                if len (meancc75)==0:
                    meancc75.append(0)
                if len(fragfitheli)==0:
                    fragfitheli.append(0)
                if len(fragfitother)==0:
                    fragfitother.append(0)
                if len(sequenceidentity)==0:
                    sequenceidentity.append(0)
                if len(seqidentfull)==0:
                    seqidentfull.append(0)
                if len(fragfitbeta)==0:
                    fragfitbeta.append(0)
                if len(fragfithelior)==0:
                    fragfithelior.append(0)
                if len(fragfitbetaor)==0:
                    fragfitbetaor.append(0)
                if len(fragitotheror)==0:
                    fragitotheror.append(0)
                if len(linkitonlybetter)==0:
                    linkitonlybetter.append(0)
                if len(fragfitonlybetter)==0:
                    fragfitonlybetter.append(0)
                if len(fragfitfirst)==0:
                    fragfitfirst.append(0)
                if len(linkitfirst)==0:
                    linkitfirst.append(0)
                if len(linkitfirstbad)==0:
                    linkitfirstbad.append(0)    
                if  len(fragfitfirstbad)==0:
                    fragfitfirstbad.append(0)
                if len(linkitfirstbadf)==0:
                    linkitfirstbadf.append(0)
                if len(fragfitfirstbadf)==0:
                    fragfitfirstbadf.append(0)
                    
                if len(helixlinkit)==0:
                    helixlinkit.append(0)
                if len(helixfragfit)==0:
                    helixfragfit.append(0)
                if len(sheetlinkit)==0:
                    sheetlinkit.append(0)
                if len(sheetfragfit)==0:
                    sheetfragfit.append(0)
                if len(fullsheetlinkit)==0:
                    fullsheetlinkit.append(0)
                if len(fullsheetfragfit)==0:
                    fullsheetfragfit.append(0)
                if len (tmscounter)==0:
                    tmscounter.append(0)
                if len (tmscounterho)==0:
                    tmscounterho.append(0)
                if len (seqidfffull)==0:
                    seqidfffull.append(0)
                if len (fsnohomoid)==0:
                    fsnohomoid.append(0)
                if len(seqidoptimal)==0:
                    seqidoptimal.append(0)
                if len(cc75)==0:
                    cc75.append(0)
                with open (t_test, 'w') as tfile:
                    st,sp = scipy.stats.ttest_rel(wodensnh,bestrmsdsnh)
                    t_value="%s:%s%s" % ("t-value",st,'\n')
                    p_value= "%s:%s" % ('p-value',sp)
                    tfile.write (t_value)
                    tfile.write(p_value)
                    print l,st,sp
                    xyz=0
                    if sp <0.05:
                        print 'yeeha'
                        xyz+=1
                    
            
                #print fragfitonlybetter,linkitonlybetter
                #break
               # print 'looplaenge:' ,l
               # print 'summe ordiff75:', sum(ordiff75)
               # print 'length ordiff75:', len(ordiff75)
               # print 'geteilt:',sum(ordiff75)/len(ordiff75)
                loopline= "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % (l,sum(bestrmsds)/len(bestrmsds),np.std(bestrmsds),
                                                         sum(wodens)/len(wodens),np.std(wodens),
                                                         sum(bestrmsdsnh)/len(bestrmsdsnh),np.std(bestrmsdsnh),
                                                         sum(wodensnh)/len(wodensnh),np.std(wodensnh),
                                                         sum(bestrmsdnhli)/len(bestrmsdnhli),sum(fragfithomo90)/len(fragfithomo90),np.std(bestrmsdnhli),
                                                         sum(wodensnhli)/len(wodensnhli),np.std(wodensnhli),
                                                         len(wodensnh),sum (optimalrmsd)/len(optimalrmsd),sum(meancorrel)/len(meancorrel),
                                                         sum(extragut)/len(extragut),sum(bestfiveall)/len(bestfiveall),
                                                         sum(besttenall)/len(besttenall),sum(linkitheli)/len(linkitheli),
                                                         sum(linkitother)/len(linkitother),sum(rank)/len(rank),
                                                         sorted(bestrmsdnhli)[len(bestrmsdnhli)//2],sum(tmsnhli)/len(tmsnhli),
                                                         sum(result50)/len(result50),sum(result20)/len(result20),sum(ori_meancorrel)/len(ori_meancorrel),
                                                         sum(or_diffcc)/len(or_diffcc),sum(ordiff75)/len(ordiff75),sum(ordiffh75full)/len(ordiffh75full),sum(diffbestfive)/len(diffbestfive),
                                                         sum(meancc75)/len(meancc75),sum(fragfitheli)/len(fragfitheli),sum(fragfitother)/len(fragfitother),sum(fragfitbeta)/len(fragfitbeta),
                                                         sum(sequenceidentity)/len(sequenceidentity),sum(seqidentfull)/len(seqidentfull),sum(fragfithelior)/len(fragfithelior),sum(fragfitbetaor)/len(fragfitbetaor),
                                                         sum(fragitotheror)/len(fragitotheror),sucounter/len(bestrmsds),sucounternh/len(bestrmsds),sum(linkitonlybetter)/len(linkitonlybetter),
                                                         sum(fragfitonlybetter)/len(fragfitonlybetter),sum(linkitfirst)/len(linkitfirst),sum(fragfitfirst)/len(fragfitfirst),
                                                         sum (linkitfirstbad)/len(linkitfirstbad),sum(fragfitfirstbad)/len(fragfitfirstbad),sum (linkitfirstbadf)/len(linkitfirstbadf),sum(fragfitfirstbadf)/len(fragfitfirstbadf),badf,
                                                         sum(helixlinkit)/len(helixlinkit),sum(helixfragfit)/len(helixfragfit),len(helixfragfit),sum(sheetlinkit)/len(sheetlinkit),sum(sheetfragfit)/len(sheetfragfit),len(sheetfragfit),
                                                         sum(fullsheetlinkit)/len(fullsheetlinkit),sum(fullsheetfragfit)/len(fullsheetfragfit),len(fullsheetfragfit),sum(tmscounter)/len(tmscounter),sum(tmscounterho)/len(tmscounterho),sum(seqidfffull)/len(seqidfffull),sum(fsnohomoid)/len(fsnohomoid),sum(seqidoptimal)/len(seqidoptimal),sum(cc75)/len(cc75),'\n')
                #print np.std(bestrmsdsnh) 
                analyse.write(loopline)
                diffile.write('\n')
                valu.write(str(l))
                valu.write(',')
                valu.write('\n')
                for j in (bestrmsdnhli):
                    valu.write(str(j))
                    valu.write(',')
                    valu.write('\n')
        print xyz
class Analyse_DSSP (PyTool):
    args= [
        _( "pdb_file", type="file" )
        ]
    #out = [
    #    _( "dssp_file", file="out.dssp" )
    #    ]
    def func (self , *args, **kwargs ):
        Dssp(self.pdb_file)
        dssp_name='%s%s' % (self.pdb_file[:-3],'dssp')
        #print dssp_name
        l=0
        h=0
        b=0
        gs=0
        x=parse_dssp(dssp_name)
        for key in x[1:]:
            gs+=1
            if type(key) == None:
                print 'huray'
                x[key].delete
            try:
                #print key[3]
                if key[3]==' ':
                    l+=1
                if key[3] in ['G','H','I']:
                   h+=1
                if key[3] in ['B','E','S','T']:
                    b+=1
            except:
                continue
            #continue
        print "loop_residues", (l/gs)*100 , 'gesamt', gs, "heli", (h/gs)*100, 'sheet', (b/gs)*100, b+h+l
        return x
class SpiderAnalyse3(PyTool):
   args = [
   _( "dataset_dir", type="dir" ),
   _( "loop_length", type="int" ),
   _( "homologue_list", type="file" ),
   _("analysis_flag", type="int")
   ]
   out = [
    _( "stat_file", file="statistics.csv" ),
    _( "res_file", file="overall.txt" )
   ]
   def func (self):

       with open ('loopanaresu.csv','w') as analyse, open ('values2.csv', 'w') as valu:
           title="%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % ("looplength","fragfit","std rmsd with density",
                                                                  "rmsd linkit","std rmsd without density",
                                                                  "rmsd fragfit no homologues","std fragfit no homologues",
                                                                  "rmsd linkit no homologs","std rmsd linkit no homologs",
                                                                  "rmsd fragfit no linker homologues","std fragfit no linker homologues",
                                                                  "rmsd linkit no linker homologs","std linkit no linker homologs",
                                                                  "Anzahl Suchen","optimal","mean_correl","extragut","beste five","linkit_helical","linkit_other","mean rank",'\n')
           analyse.write(title)
           correlation={
           '5':0.568783917098,
           '6':0.56042945098,
           '7':0.523105761658,
           '8':0.463180293532,
           '9':0.41430411,
           '10':0.367170015707,
           '11':0.329382698022,
           '12':0.298228513661,
           '13':0.250647477612,
           '14':0.222953401736,
           '15':0.199863570246,
           '16':0.174128295472,
           '17':0.167516662037,
           '18':0.10233948813,
           '19':0.0990273932505,
           '20':0.0832893506279,
           '21':0.0721325818392,
           '22':0.0700303510774,
           '23':0.0653087346573,
           '24':0.0603087346573,
           '25':0.05819717193,
           '26':0.0612093805726,
           '27':0.0621014396482,
           '28':0.0568244603245,
           '29':0.0757444702719,
           '30':0.0770906372842,
           '31':0.0669438991446,
           '32':0.0775876085589,
           '33':0.0793356798089,
           '34':0.0514788622987,
           '35':0.0514788622987
           
           }
           for l in range(5,35,1):
               scorr=correlation[str(l)] #standartdcorreltion fuer die looplaenge
               corrband1=scorr+scorr/5
               corrband2=scorr-scorr/5
               print corrband2, corrband1
               #break
               #corrrange=range(corrband2,corrband1)
               chainfolders=['2','3']
               bestrmsds=[]
               bestrmsdsnh=[] #bestrmsdsno homologues
               wodens=[]
               wodensnh=[] #without density no homologues
               wodensnhli=[]
               bestrmsdnhli=[] #no homologue linker
               meancorrel=[]    
               optimalrmsd=[]
               extragut=[]
               bestfiveall=[]
               linkitheli=[]
               linkitother=[]
               rank=[]
               hc=0
               hc2=0
               
               for g in chainfolders:
                   print g
                   loopfolder=os.path.join(self.dataset_dir,g,str(l))
                   if os.path.isdir(loopfolder):
                       fi=sorted((os.listdir(loopfolder)), key= lambda x: int(x.split('_')[-1]))
                       #folders=os.listdir(self.dataset_dir)
                       z=2
                       homo = [line.rstrip() for line in open(self.homologue_list)]
               
                       with open ("analysis.csv",'w') as ana:
                           
                           print fi
                           for x in fi:
                               y= int(filter(lambda x: x.isdigit(),x[-4:]))
                               if y<999:
                                   name="%s.%s" % (x,'pdb')
                                   ori=os.path.join(loopfolder,x,name)
                                   orinpdb=numpdb.NumPdb(ori, {
                                                   "phi_psi": False,
                                                   "sstruc": False,
                                                   "backbone_only": True,
                                                   "protein_only": False,
                                                   "detect_incomplete": False,
                                                   "configuration": False,
                                                   "info": False
                                                   })
                                   result={}
                                   resultnh={}
                                   #resultnh={}
                                   resultlist=[]
                                   seq1=orinpdb.sequence()
                                   loops=os.path.join(loopfolder,x,'loop_correl/crosscorrelation/loops')
                                   ccsort=os.path.join(loopfolder,x,'loop_correl/crosscorrelation/ccsort.cpv')
                                   linkertxt=os.path.join(loopfolder,x,'link_it/edited_linker.txt')
                                   if os.path.isfile(ccsort):
                                       with open (ccsort,'r') as hurz:
                                           lines=hurz.readlines()
                                           nr=lines[1].split()[0]
                                           corr=lines[1].split()[2]
                                           meancorrel.append(float(corr))
                                           #hurz=','.join(lines[1].split())
                                           #print nr
                                          # break
                                       with open (linkertxt, 'r') as linki:
                                           lineslinker=linki.readlines()
                                       
                                       h=os.listdir(loops)
                                       for i in h:
                                           if i.endswith('.pdb'):
                                               #print i
                                               
                                               lpath=os.path.join(loops,i)
                                               npdb = numpdb.NumPdb(lpath,{
                                                           "phi_psi": False,
                                                           "sstruc": False,
                                                           "backbone_only": True,
                                                           "protein_only": False,
                                                           "detect_incomplete": False,
                                                           "configuration": False,
                                                           "info": False
                                                           })
                                               #print orinpdb['xyz']
                                              # print lpath
                                               lnr=int(i[0:3]) *4+1
                                               seqnr=int(i[0:3]) *4
                                               #print orinpdb['xyz']
                                               #print npdb['xyz']
                                               #Dssp(lpath)
                                               rootm2=rmsd(orinpdb['xyz'],npdb['xyz'])
                                               seq2=lineslinker[seqnr].strip()
                                               
                                               #print seq1
                                               #sequence identity
                                               si=0
                                               for sf in range(0, len(seq1),1):
                                                   if seq1[sf]==seq2[sf+1]:
                                                       si+=1
                                               sqi=(si/len(seq1))*100
                                               #print seq1,"and",seq2[1:-1]
                                               #print sqi
                                                   
                                               
                                                   
                                               #print rootm2
                                               result[i]=[rootm2]
                                               result[i].append(lineslinker[lnr].strip())
                                               result[i].append(sqi)
                                               resultlist.append(rootm2)
                                               
                                           else:
                                               continue
                                           
                                       #den result dictonary auf eintraege resuzieren die nicht auf der Liste stehen:
                                       resultnh=result.copy()
                                       for key in result:
                                           if result[key][1].upper() in homo:
                                               del resultnh[key]
                                       #einen dictonary erstellen wo die homologen linker nicht drin sind:
                                       resultnhli=result.copy()
                                       for key in result:
                                           if result[key][2] >=75:
                                               del resultnhli[key]
                                      # print "ganzer dict", len(result)
                                      # print "homo liste", len(resultnh)
                                      # print "homo linker", len(resultnhli)
                                       pdbkey="%s_%s.%s" % (nr.zfill(3),'bb','pdb')
               
                                       
                                       #abchecken ob 
                                       #if result[pdbkey][1].upper() in homo:
                                       #    hc+=1
                                       #    print 'yes'
                                       #else:
                                       #    continue
               
                                       best=sorted(resultlist)
                                       #print result
                                      # print "\n",x,"\n"
                                       #rmsd vonn loop nr 1
                                      # print "without desnity:" ,result['001_bb.pdb']
                                       #rmsd vom hoechst geranketem loop
                                      # print "best:",pdbkey, result[pdbkey]
                                       #der optimale rmsd: also einfach der kleinste aller 100
                                       #print "optimal",min(result, key=result.get),best[0]
                                       nr_opt=min(result, key=result.get)
                                       
                                       i=1
                                       
                                       # an welcher Stelle is loop nr 1 gerankt
                                       
                                       for c,line in enumerate (lines, 1):
                                           gu=line.split()
                                           if gu[0]==str(nr_opt)[1:3]:
                                               z=c
                                               rank.append(c-1)
                                              # print "was ranked on number:", c-1
                                           i+=1
                                       #der schlechstest moegliche rmsd, also der letzte in der sortierten liste
                                      # print "worst",best[-1]
                                       if pdbkey==min(result, key=result.get):
                                           opti=1
                                       else:
                                           opti=0
                                       if result[pdbkey][0]<result['001_bb.pdb'][0]:
                                           better=1
                                       else:
                                           better=0
                                       rline="%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % (x,result['001_bb.pdb'][0],pdbkey,result[pdbkey][0],min(result, key=result.get),best[0],z-1,best[-1],opti,better,'\n')
                                       bestrmsds.append(result[pdbkey][0])
                                       wodens.append(result['001_bb.pdb'][0])
                                       optimalrmsd.append(best[0])
                                       ###nur extragute
                                       if corrband2<float(corr) <corrband1 and result[pdbkey][2]<75:
                                           extragut.append(result[pdbkey][0])
                                       #
                                       ana.write(rline)
                                       
                                       
                                       ##############
                                       
                                       #schleife ueber die loops
                                       
                                       for y in range (1,100,1):
                                           loop="%s_%s.%s" % (str(y).zfill(3),'bb','pdb')
                                          # print loop
                                           if loop in resultnh:
                                               wodensnh.append(result[loop][0])
                                               break
                                           else:
                                               continue
                                           
                                       # das ganze fuer mit density
                                       
                                       for line in lines[1:]:
                                           #lines.next()
                                           no=line.split()[0]
                                           lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                           
                                          # print 'richtig so?', result[lonam][1].upper()
                                           if result[lonam][1].upper() not in homo:
                                               #print 'bloo'
                                               bestrmsdsnh.append(result[lonam][0])
                                               break
                                           else:
                                               continue
                                       #homologe linker aussortieren:
                                       for y in range (1,100,1):
                                           loop="%s_%s.%s" % (str(y).zfill(3),'bb','pdb')
                                          # print loop
                                           if loop in resultnhli:
                                               wodensnhli.append(result[loop][0])
                                               hfj=os.path.join(loops,loop)
                                               hjd='%s%s' % (hfj[:-3],'dssp')
                                               Dssp(hfj,output_dir=loops)
                                               dsrec=parse_dssp(hjd)
                                             #  print hjd
                                               
                                               #with open ('hjd') as ds:
                                               #    gd=ds.readlines()
                                               #    for line in gd:
                                               #        if line.split() []
                                               if any('H'==e[3]for e in dsrec[1:])==True:
                                                   print 'found helical loop'
                                                   linkitheli.append(result[loop][0])
                                               else:
                                                   linkitother.append(result[loop][0])
                                               #for i in dsrec[1:]:
                                               #    if i[3]=='H':
                                               #        print i[3]
                                               #print "dssp", dsrec[1]
                                               break
                                           else:
                                               continue
                                           
                                       
                                       for line in lines[1:]:
                                           #lines.next()
                                           no=line.split()[0]
                                           lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                           if result[lonam][2]<75:
                                               bestrmsdnhli.append(result[lonam][0])
                                               break
                                           else:
                                               continue
                                       
                                       #die besten five
                                       bestfive=[]
                                       for line in lines[1:]:
                                           #lines.next()
                                           no=line.split()[0]
                                           lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                           if result[lonam][2]<75:
                                               bestfive.append(float(result[lonam][0]))
                                               if len(bestfive)==5:
                                                   #print sorted(bestfive)
                                                   bestfiveall.append(sorted(bestfive)[0])
                                                   break
                                           else:
                                               continue
                                       print '.'
               
               if len(linkitheli)==0:
                   linkitheli.append(0)
               if len(extragut)==0:
                   extragut.append(0)
               if len(bestrmsds)==0:
                   bestrmsds.append(0)
               if  len(wodens)==0:
                   wodens.append(0)
               if len(optimalrmsd)==0:
                   optimalrmsd.append(0)
               if len(bestrmsdsnh)==0:
                   bestrmsdsnh.append(0)
               if len(wodensnh)==0:
                   wodensnh.append(0)
               if len(bestrmsdnhli)==0:
                   bestrmsdnhli.append(0)
               if len(wodensnhli)==0:
                   wodensnhli.append(0)
               if  len(meancorrel)==0:
                   meancorrel.append(0)
               if len(bestfiveall)==0:
                   bestfiveall.append(0)
               if len(linkitother)==0:
                   linkitother.append(0)
               if len(rank)==0:
                   rank.append(0)
               print "besten five",bestfiveall                     
               print hc
               print 'mean',sum(bestrmsds)/len(bestrmsds),'anzahl:',len(bestrmsds)
               print 'mean wo dens', sum(wodens)/len(wodens),'anzahl:',len(wodens)
               print 'best possible', sum (optimalrmsd)/len(optimalrmsd)
               print 'mean no homologues',sum(bestrmsdsnh)/len(bestrmsdsnh),'anzahl:',len(bestrmsdsnh)
               print 'mean wo dens no homologue', sum(wodensnh)/len(wodensnh),'anzahl:',len(wodensnh)
               print 'extragut', sum(extragut), 'anzahl', len(extragut)
               loopline= "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % (l,sum(bestrmsds)/len(bestrmsds),np.std(bestrmsds),
                                                        sum(wodens)/len(wodens),np.std(wodens),
                                                        sum(bestrmsdsnh)/len(bestrmsdsnh),np.std(bestrmsdsnh),
                                                        sum(wodensnh)/len(wodensnh),np.std(wodensnh),
                                                        sum(bestrmsdnhli)/len(bestrmsdnhli),np.std(bestrmsdnhli),
                                                        sum(wodensnhli)/len(wodensnhli),np.std(wodensnhli),
                                                        len(wodensnh),sum (optimalrmsd)/len(optimalrmsd),sum(meancorrel)/len(meancorrel),sum(extragut)/len(extragut),sum(bestfiveall)/len(bestfiveall),sum(linkitheli)/len(linkitheli),sum(linkitother)/len(linkitother),sum(rank)/len(rank),sucounter/len(bestrmsds),'\n')
               print np.std(bestrmsdsnh) 
               analyse.write(loopline)
               valu.write(str(l))
               valu.write(',')
               valu.write('\n')
               for j in (bestrmsdnhli):
                   valu.write(str(j))
                   valu.write(',')
                   valu.write('\n')
                   
class LinkItRMSD(PyTool):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "dataset_dir", type="dir" ),
        _( "res1", type="sele" ),
        _( "res2", type="sele" )
           ]
    
    def func( self):
        with open ('rmsd.csv','w') as rms:
            print 'hallo'
            res1=self.res1['resno']
            res2=self.res2['resno']
            cha=self.res1['chain']
            npdb=numpdb.NumPdb( self.pdb_file,{"backbone_only": True})
            oripdb=npdb.sele(resno=[res1,res2],chain=cha, invert=False)
            orinpdb=npdb.copy(sele=oripdb)
            npdb.write('ori.pdb',sele=oripdb)
            #loopdir=os.path.join()
            z=os.path.join(self.dataset_dir,'loops_repair')
            h=os.listdir(z)
            #print h
            for key in orinpdb:
                print key
                print key['atomname']
            
            hurz=[h for h in os.listdir(z) if any ([h.endswith('.pdb')]) ]
            fi=sorted(hurz, key= lambda x: int(x.split('_')[0]))
            print fi
            for i in fi:
                print i
               # if i.endswith('.pdb'):
                x=os.path.join(self.dataset_dir,'loops_repair',i)
                loopnpdb=numpdb.NumPdb(x)
                rootm2=rmsd(orinpdb['xyz'],loopnpdb['xyz'])
                line="%s,%s%s" % (i,rootm2,'\n')
                rms.write(line)
                print rootm2
class RepairLoops(PyTool):
    args = [ 
        _( "dataset_dir", type="dir" )
        ]
    out = [
        _( "loop_dir", dir="loops_repair")]
    def func(self):
        for l in range(5,35,1):
            chainfolders=os.listdir(self.dataset_dir)
            for g in chainfolders:
                loopfolder=os.path.join(self.dataset_dir,g,str(l))
                if os.path.isdir(loopfolder):
                    loopfolders=os.listdir(loopfolder)
                    
                    for tu in loopfolders:
                        loopcorrel=os.path.join(loopfolder,tu,'loop_correl')
                        if os.path.exists(loopcorrel):
                            z=os.path.join(loopfolder,tu,'loop_correl/crosscorrelation/loops')
                            zi=os.path.join(loopfolder,tu,'loop_correl/crosscorrelation/loops_repair')
                            if not os.path.exists(zi):
                                os.mkdir(zi)
                                hurz=[h for h in os.listdir(z) if any ([h.endswith('.pdb')]) ]
                                for i in hurz:
                                    x=os.path.join(z,i)
                                    y=os.path.join(zi,i)
                                    print y
                                    with open (x,'r')as huhu,open( y,'w') as buhu:
                                        gu=huhu.readlines()
                                        for i in range (0,len(gu)-2,4):
                                            atom1="%7s" % str(i+1)
                                            line1=gu[i][0:4] + atom1 +gu[i][11:]
                                            atom2="%7s" % str(i+2)
                                            line2=gu[i+2][0:4] + atom2 +gu[i+2][11:]
                                            atom3="%7s" % str(i+3)
                                            line3=gu[i+1][0:4] + atom3 +gu[i+1][11:]
                                            atom4="%7s" % str(i+4)
                                            line4=gu[i+3][0:4] + atom4 +gu[i+3][11:]
                                            print line4
                                            buhu.write(line1)
                                            buhu.write(line2)
                                            buhu.write(line3)
                                            buhu.write(line4)
                                        buhu.write('END')
                                        buhu.close()
                        else:
                            continue
class RepairLinkItLoops(PyTool):
    args = [ 
    _( "dataset_dir", type="dir" )
    ]
    out = [
        _( "loop_dir", dir="loops_repair")]
    def func(self):
        z=os.path.join(self.dataset_dir,'loops')
        zi=os.path.join(self.dataset_dir,'loops_repair')
        if not os.path.exists(zi):
                os.mkdir(zi)
                hurz=[h for h in os.listdir(z) if any ([h.endswith('.pdb')]) ]
                for i in hurz:
                    x=os.path.join(z,i)
                    y=os.path.join(zi,i)
                    print y
                    with open (x,'r')as huhu,open( y,'w') as buhu:
                        gu=huhu.readlines()
                        for i in range (0,len(gu)-2,4):
                            atom1="%7s" % str(i+1)
                            line1=gu[i][0:4] + atom1 +gu[i][11:]
                            atom2="%7s" % str(i+2)
                            line2=gu[i+2][0:4] + atom2 +gu[i+2][11:]
                            atom3="%7s" % str(i+3)
                            line3=gu[i+1][0:4] + atom3 +gu[i+1][11:]
                            atom4="%7s" % str(i+4)
                            line4=gu[i+3][0:4] + atom4 +gu[i+3][11:]
                            print line4
                            buhu.write(line1)
                            buhu.write(line2)
                            buhu.write(line3)
                            buhu.write(line4)
                        buhu.write('END')
                        buhu.close()
class RepairLoops2(PyTool):
    args = [ 
        _( "dataset_dir", type="dir" )
        ]
   #out = [
    #    _( "loop_dir", dir="loops_repair")]
    def func(self):
        chainfolders=os.listdir(self.dataset_dir)
        for g in chainfolders:
        #    loopfolder=os.path.join(self.dataset_dir,g,str(l))
        #    if os.path.isdir(loopfolder):
        #        loopfolders=os.listdir(loopfolder)
        #        
        #        for tu in loopfolders:
            folder=os.path.join(self.dataset_dir)
        #            if os.path.exists(loopcorrel):
        
            z=os.path.join(folder,g,'loops')
            zi=os.path.join(folder,g,'loops_repair')
            if not os.path.exists(zi):
                os.mkdir(zi)
                hurz=[h for h in os.listdir(z) if any ([h.endswith('.pdb')]) ]
                for i in hurz:
                    x=os.path.join(z,i)
                    y=os.path.join(zi,i)
                    print y
                    with open (x,'r')as huhu,open( y,'w') as buhu:
                        gu=huhu.readlines()
                        for i in range (0,len(gu)-2,4):
                            atom1="%7s" % str(i+1)
                            line1=gu[i][0:4] + atom1 +gu[i][11:]
                            atom2="%7s" % str(i+2)
                            line2=gu[i+2][0:4] + atom2 +gu[i+2][11:]
                            atom3="%7s" % str(i+3)
                            line3=gu[i+1][0:4] + atom3 +gu[i+1][11:]
                            atom4="%7s" % str(i+4)
                            line4=gu[i+3][0:4] + atom4 +gu[i+3][11:]
                            print line4
                            buhu.write(line1)
                            buhu.write(line2)
                            buhu.write(line3)
                            buhu.write(line4)
                        buhu.write('END')
                        buhu.close()
            else:
             continue
class SuperLooperAnalyse (PyTool):
    args = [
        #_( "pdb_file", type="file", ext="pdb" ),
        _( "linklist", type="file", ext="txt" ),
        _( "dataset_dir", type="dir" )
        #_( "res1", type="sele" ),
       # _( "res2", type="sele" )
           ]
        
    def func(self):
        result={}
        chainfolders=os.listdir(self.dataset_dir)
        with open ('analysis.cpv','w') as ana,open('analysis_good.csv','w')as anago:
            for g in chainfolders:
                fname="%s_%s" %  (g[0:4],'linker.txt')
                pdb_name="%s.%s" % (g[0:4],'pdb')
                linkertxt=os.path.join(self.dataset_dir,g,fname)
                print fname
                linkitfile=os.path.join(self.dataset_dir,g,'linkit.json')
                with open(self.linklist,'r') as linki,open('fails.txt','w') as fail:
                    gul=linki.readlines()
                    for line in gul[1:-1]:
                        pdb=line.split( )[0]
                        as1=line.split( )[1]
                        as2=line.split( )[2]
                        u="%s_%s_%s" % (pdb,as1,as2)
                        result[u]=[line.split( )[11]]
                        result[u].append(line.split( )[5])
                        result[u].append(line.split( )[13])
                        result[u].append(line.split( )[4])
                        seq=line.split( )[3]
                with open (linkertxt, 'r') as linki, open (linkitfile,'r')as jasi:
                    lineslinker=linki.readlines()
                    stems=json.load(jasi)
                    stem1=stems[0] [1]
                    stem2=stems[0] [2]
                    res1=stem1.split(':')[0]
                    res2=stem2.split(':')[0]
                    chain=stem2.split(':')[1]
                    pdb_file=os.path.join('downloadpdb',pdb_name)
                    npdbori=numpdb.NumPdb(pdb_file,
                                          {"phi_psi": False,
                                "sstruc": False,
                                "backbone_only": True,
                                "protein_only": False,
                                "detect_incomplete": False,
                                "configuration": False,
                                "info": False
                            })
                    test=npdbori.copy(resno=[int(res1)+1,int(res2)-1],chain=chain)
                    print res1,res2,chain
                    
                    loops=os.path.join(self.dataset_dir,g,'loops_repair')
                    h=os.listdir(loops)
                    fi=sorted(h, key= lambda x: int(x.split('_')[0]))
                    seq=test.sequence()
                    for i in fi:
                        if i.endswith('.pdb'):
                ##print i
                            lino=int(i[0:3])
                            lpath=os.path.join(loops,i)
                            npdb = numpdb.NumPdb(lpath,{
                                "phi_psi": False,
                                "sstruc": False,
                                "backbone_only": True,
                                
                                "protein_only": False,
                                "detect_incomplete": False,
                                "configuration": False,
                                "info": False
                            })
                            npdbca=npdb.copy()
                            
                            lnr=int(i[0:3]) *4+1
                            seqnr=int(i[0:3]) *4
                            
                            pdbori=lineslinker[lnr].strip()
                            seq2=lineslinker[seqnr].strip()
                            print pdbori, seq,seq2[1:-1]
                           # print seqnr
                            out="%s_%s.%s" % (pdb_name[:-4],'orilinker','pdb')
                            outf=os.path.join(self.dataset_dir,g,out)
                            test.write(outf)
                            backbone = ( ' N  ',' C  ', ' CA ',' O  ' )
                            sele1={"atomname":backbone}
                            sele2={"atomname":backbone}
                            l=npdb.sele(atomname=backbone)
                           #a,b=numpdb.superpose(test,npdb,sele1,sele2)
                            
                            
                            #rootm5=rmsd[a['xyz'],b['xyz']]
                            #print 'superposed', rootm5
                            rootm2=rmsd(test['xyz'],npdbca['xyz'])
                            if rootm2==0:
                                #print 'yeha'
                                #print i
                                continue
                            else:
                                sqi=0
                                si=0
                                for sf in range(0, len(seq),1):
                                    if seq[sf]==seq2[sf+1]:
                                        si+=1
                                sqi=(si/len(seq))*100
                                backbone = ( ' N  ',' C  ', ' CA ',' O  ' )
                                sele1={"atomname":backbone}
                                sele2={"atomname":backbone}
                                l=npdb.sele(atomname=backbone)
                                a,b=numpdb.superpose(test,npdb,sele1,sele2)
                                rootm3=rmsd(test['xyz'],npdbca['xyz'])
                                print i,rootm2,rootm3#sqi, result[g][0],result[g][1]
                                line_ana="%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % (g,i[0:3],pdbori,round(rootm2,2),result[g][0],result[g][1],round(sqi,2),len(seq2),result[g][2],rootm3,'\n')
                                #print line_ana
                                ana.write(line_ana)
                                if int(result[g][3])==2:
                                   anago.write (line_ana)
                                else:
                                    continue
                                break
class ShortenLoops (PyTool):
    args = [
        #_( "pdb_file", type="file", ext="pdb" ),
        _( "dataset_dir", type="dir" )
        #_( "res1", type="sele" ),
       # _( "res2", type="sele" )
           ]
        
    def func(self):
        chainfolders=os.listdir(self.dataset_dir)
        for g in chainfolders:
            z=os.path.join(self.dataset_dir,g,'loops')
            zi=os.path.join(self.dataset_dir,g,'loops_repair')
            hurz=[h for h in os.listdir(z) if any ([h.endswith('.pdb')]) ]
            for i in hurz:
                x=os.path.join(z,i)
                #y=os.path.join(zi,i)
                npdb=numpdb.NumPdb(x)
                seq=npdb.sequence()
                l=len(seq)
                res1=npdb.get('resno')[0]
                res2=res1+l
                print res1, res2
                test=npdb.copy(resno=[res1+1,res2-1])
                test.write(x)
                
            hurz=[h for h in os.listdir(zi) if any ([h.endswith('.pdb')]) ]
            for i in hurz:
                x=os.path.join(zi,i)
                #y=os.path.join(zi,i)
                npdb=numpdb.NumPdb(x)
                seq=npdb.sequence()
                l=len(seq)
                res1=npdb.get('resno')[0]
                res2=res1+l
                print res1, res2
                test=npdb.copy(resno=[res1+1,res2-1])
                test.write(x)
class MultiLooporicc( PyTool ):
    args = [
        _( "boxsize",type="float"),
        _( "boxsize2",type="float"),
        _( "boxsize3",type="float"),
        _("pdb_file", type="file"),
        _( "dataset_dir", type="dir" )
        ]
    def func ( self, *args, **kwargs ):
            #l=25
        for l in range(5,36,1):
            #l=19
            chainfolders=os.listdir(self.dataset_dir)
            for g in chainfolders:
                
                loopfolder=os.path.join(self.dataset_dir,g,str(l))
                if os.path.isdir(loopfolder):
                   
                    
                    print 'G',g
                    fi=sorted((os.listdir(loopfolder)), key= lambda x: int(x.split('_')[-1]))#code
                    for x in fi:
                        name="%s.%s" % (x,'pdb')
                        oriccfile=os.path.join(loopfolder,x,'oricrosscorrelation.cpv')
                        if os.path.exists(oriccfile):
                            continue
                        else:
                            res1=x.split('_')[0]
                            res2=x.split('_')[1]
                            
                            npdb=numpdb.NumPdb( self.pdb_file)
                            
                            test=npdb.copy(resno=[int(res1),int(res2)],chain=g)
                            out=os.path.join(loopfolder,x)
                            namefull= "%s/%s_%s.%s" % (out,res1,res2,'pdb')
                            test.write(namefull)
                            subs=[h for h in os.listdir(out) if any ([h.startswith('loop_correl')]) ]
                            for fo in subs:
                                #print res1,res2,fo
                                ori=os.path.join(loopfolder,x,name)
                                mapori=os.path.join(loopfolder,x,fo,'delete_filled_densities','usermap.cpv')
                                mask=os.path.join(loopfolder,x,fo,'crosscorrelation','mask.cpv')
                                #print maske
                                #print ori
                                shiftfile=os.path.join(loopfolder,x,fo,'shift','shift.cpv')
                                eckenfile=os.path.join(loopfolder,x,fo,'box','ergebnisse.cpv')
                                b_out=os.path.join(loopfolder,x,name)
                                
                                outd=os.path.join(loopfolder,x,fo)
                                #print outd
                                ccsortf=os.path.join(loopfolder,x,fo,'crosscorrelation','ccsort.cpv')
                                #pdb_sh=os.path.join(loopfolder,x,fo,'edited.pdb')
                                if os.path.isfile(eckenfile):
                                    print 'eckenfile'
                                    if os.path.isfile(ccsortf):
                                        print 'calculate',loopfolder,x,fo,name
                                        print outd
                                        with open (eckenfile, 'r') as ecki,open (shiftfile,'r') as shifti:
                                            guja=ecki.readlines()
                                            h=guja[1].split()
                                            obenlinks1=h[3]
                                            obenlinks2=h[4]
                                            obenlinks4=h[6]
                                            ps=h[7]
                                            reso=h[8]
                                            obs=h[2]
                                            #print h[4]
                                           # print out
                                            #print
                                            blue=shifti.readlines()
                                            f=blue[1].split()
                                            print f[2]
                                            sh1=float(f[2])
                                            sh2=float(f[3])
                                            sh3=float(f[4])
                                            #PdbEdit(ori,shift=[sh1,sh2,sh3],output_dir=outd)
                                            
                                        Calcoricc(self.boxsize,self.boxsize2,self.boxsize3,obs,mapori,mask,ori,obenlinks1,obenlinks2,obenlinks4,reso,ps,output_dir=outd)
                                        #os.remove(b_out)
                                        #os.rename(pdb_sh,b_out)
class MultiSpiderAnalyse2(PyTool):
    args = [
    _( "dataset_dir", type="dir" ),
    _( "loop_length", type="int" ),
    _( "homologue_list", type="file" ),
    _( "header_pdb", type="file" ),
    _("analysis_flag", type="int")
    ]
    out = [
     _( "stat_file", file="statistics.csv" ),
     _( "res_file", file="overall.txt" )
    ]
    def func (self):
        print "HHAAAAAAAAAAAAAAAALLLLLLLLLLOOOOOOOOOOOO"
        for cor in range(1,7,1):
            #print 'runde', cor
            helix={}
            sheet={}
            fullsheet={}
            sheetnames=[]
            with open (self.header_pdb) as header:
                headerlines=header.readlines()
                xyz=0
                for lines in headerlines:
                    if lines.startswith  ("HELIX"):
                        ident="%s_%s_%s" % ('HELIX',lines[19:20],lines[8:10])
                        helix[ident]=[int(lines[72:76])]
                        helix[ident].append(lines[8:10])
                        helix[ident].append(lines[31:32])
                        helix[ident].append(int(lines[22:25]))
                        helix[ident].append(lines[34:37])
                    if lines.startswith  ("SHEET"):
                        #print int(lines[34:37])
                        xyz+=1
                        #hurzi=lines.split()
                        #print hurzi[4]
                        ident=lines[12:14],lines[8:10]#"%s_%s_%s.%s" % ("SHEET",lines[21:22],lines[8:10],xyz)
                        sheetlen=(int(lines[34:37])-int(lines[23:26]))+1
                        #print sheetlen
                        sheet[ident]=[sheetlen]#[int(lines[72:76])]
                        sheet[ident].append(lines[8:10])
                        sheet[ident].append(lines[21:22])
                        sheet[ident].append(int(lines[23:26]))
                        sheet[ident].append(int(lines[34:37]))
                        sheet[ident].append(lines[12:14])
                        sheetnames.append(lines[12:14])
            for hurz in sorted(set(sheetnames)):
                go={}
                for key in sorted(sheet.iterkeys()):
                    
                    if key[0]==hurz:
    
                        go[key]=sheet[key]
    
                for i in range (0,len(go)-1):
    
                    if go.values()[i] [2]==go.values()[i+1] [2] and abs(go.values()[i] [4]-go.values()[i+1] [4])<10 :
     
                        sheetlist=[]
                        sheetlist.extend([go.values()[i] [4],go.values()[i] [3],go.values()[i+1] [3],go.values()[i+1] [4]])
                        identfull= "%s_%s_%s" % (hurz,go.values()[i] [2],min(sheetlist))
                        fullsheet[identfull]=[go.values()[i] [2]]
                        fullsheet[identfull].append(min(sheetlist))
                        fullsheet[identfull].append(max(sheetlist))
                        fullsheet[identfull].append(max(sheetlist)-min(sheetlist))
            for key in sorted(fullsheet):
                #print key, fullsheet[key]
                continue
            ananame="%s_%s.%s" % ('g',cor,'csv')
            with  open (ananame,'w') as analyse,open ('values2.csv', 'w') as valu, open ('cases.txt', 'w') as cases, open ('bettercc.txt', 'w') as bettercc, open ('diff.csv','w') as diffile:
                title="%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % ("looplength","fragfit","std rmsd with density",
                                                                       "rmsd linkit","std rmsd without density",
                                                                       "rmsd fragfit no homologues","std fragfit no homologues",
                                                                       "rmsd linkit no homologs","std rmsd linkit no homologs",
                                                                       "rmsd fragfit no homologues 75","rmsd fragfit no homologues 90","std fragfit no linker homologues",
                                                                       "rmsd linkit no linker homologs","std linkit no linker homologs",
                                                                       "Anzahl Suchen","optimal","mean_correl","extragut","beste five","best ten","linkit_helical",
                                                                       "linkit_other","mean rank","median fragfit","tms scrore","Fragfit 50","Fragfit 20",
                                                                       'ori cc mean','Diff to oricc','diff to oricc homo75','diff 75','ccbestfive','meanccbest','fragfitheli','fragfitother','fragfitbeta','sequenceidentity no homo','sequence identity full','fragfithel ori','fragfitbeta ori','fragfitother ori','anteil perfekt','anteilperfekt nh','linkitbetter','fragfitbetter','linkitfirst','fragfitfirst',
                                                                       'linkit>5','fragfitof>5','linkitbad','fragfitbad','conuter bad',
                                                                       "helix_link_it","helix_fragfit",'anzahl helix',"sheet_linkit",'Sheet_fragfit','anzahl sheet','fullsheetlinkit','fullsheetfragfit','anzahl fullsheet','\n')
                analyse.write(title)
                rmsdval={
                '5': 1.7857142857,
                '6':2.1428571429,
                '7':2.5,
                '8':2.8571428571,
                '9':3.2142857143,
                '10':3.5714285714,
                '11':3.9285714286,
                '12':4.2857142857,
                '13':4.6428571429,
                '14':5,
                '15':5.3571428571,
                '16':5.7142857143,
                '17':6.0714285714,
                '18':6.4285714286,
                '19':6.7857142857,
                '20':7.1428571429,
                '21':7.5,
                '22':7.8571428571,
                '23':8.2142857143,
                '24':8.5714285714,
                '25':8.9285714286,
                '26':9.2857142857,
                '27':9.6428571429,
                '28':10,
                '29':10.3571428571,
                '30':10.7142857143,
                '31':11.0714285714,
                '32':11.4285714286,
                '33':11.7857142857,
                '34':12.1428571429
                  
                            }
                correlation={
                '5':0.568783917098,
                '6':0.56042945098,
                '7':0.523105761658,
                '8':0.463180293532,
                '9':0.41430411,
                '10':0.367170015707,
                '11':0.329382698022,
                '12':0.298228513661,
                '13':0.250647477612,
                '14':0.222953401736,
                '15':0.199863570246,
                '16':0.174128295472,
                '17':0.167516662037,
                '18':0.10233948813,
                '19':0.0990273932505,
                '20':0.0832893506279,
                '21':0.0721325818392,
                '22':0.0700303510774,
                '23':0.0653087346573,
                '24':0.0603087346573,
                '25':0.05819717193,
                '26':0.0612093805726,
                '27':0.0621014396482,
                '28':0.0568244603245,
                '29':0.0757444702719,
                '30':0.0770906372842,
                '31':0.0669438991446,
                '32':0.0775876085589,
                '33':0.0793356798089,
                '34':0.0514788622987,
                '35':0.0514788622987
                
                }
                for l in range(5,35,1):
                    indifile="%s_%s.%s" % (l,'file','csv')
                    #with open (indifile, 'w') as loofile:
                    #title_loo="%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % ('chain','name','taken','ori_pdb','seq_identit','rmsd','SSE flag','diffccabs','tmscore','diffreal','bestes_modell','bestrmsd','fragfitbest','fragfitbesrmsd','bestpossibleall?','bestlinkitnh','bestpossiblenh','linkitbest','fragfitbest','linkit>5','fragfitof>5','linkitbad','fragfitbad','\n')
                    #loofile.write(title_loo)
                    diffile.write(str(l) + ',')
                    scorr=correlation[str(l)] #standartdcorreltion fuer die looplaenge
                    corrband1=scorr+scorr/5
                    corrband2=scorr-scorr/5
                    #print corrband2, corrband1
                    #break
                    #corrrange=range(corrband2,corrband1)
                    chainfolders=os.listdir(self.dataset_dir)
                    bestrmsds=[]
                    bestrmsdsnh=[] #bestrmsdsno homologues
                    wodens=[]
                    wodensnh=[] #without density no homologues
                    wodensnhli=[]
                    bestrmsdnhli=[] #no homologue linker
                    tmsnhli=[]
                    meancorrel=[]
                    ori_meancorrel=[]
                    optimalrmsd=[]
                    extragut=[]
                    bestfiveall=[]
                    besttenall=[]
                    linkitheli=[]
                    linkitother=[]
                    fragfithomo90=[]
                    rank=[]
                    result50=[]
                    result20=[]
                    hc=0
                    hc2=0
                    sucounter=0
                    sucounternh=0
                    greatdiff=[]
                    or_diffcc=[]
                    ordiff75=[]
                    ordiffh75full=[]
                    diffbestfive=[]
                    meancc75=[]
                    fragfitheli=[]
                    fragfitother=[]
                    fragfitbeta=[]
                    sequenceidentity=[]
                    seqidentfull=[]
                    fragfithelior=[]
                    fragfitbetaor=[]
                    fragitotheror=[]
                    linkitonlybetter=[]
                    fragfitonlybetter=[]
                    fragfitfirst=[]
                    linkitfirst=[]
                    fragfitfirstbad=[]
                    linkitfirstbad=[]
                    fragfitfirstbadf=[]
                    linkitfirstbadf=[]
                    helixlinkit=[]
                    helixfragfit=[]
                    sheetlinkit=[]
                    sheetfragfit=[]
                    fullsheetlinkit=[]
                    fullsheetfragfit=[]
                    badf=0

                    for g in chainfolders:
                        #g='A'
                        #print g
                        loopfolder=os.path.join(self.dataset_dir,g,str(l))
                        if os.path.isdir(loopfolder):
                            fi=sorted((os.listdir(loopfolder)), key= lambda x: int(x.split('_')[-1]))
                            #folders=os.listdir(self.dataset_dir)
                            z=2
                            homo = [line.rstrip() for line in open(self.homologue_list)]
                            #print loopfolder
                            with open ("analysis.csv",'w') as ana:
                                
                                #print fi
                                for x in fi:
                                    #print cor, x
                                    y= int(filter(lambda x: x.isdigit(),x[-4:]))
                                    y1=int(filter(lambda x: x.isdigit(),x[:3]))
                                   #print y1,y
                                    if y<999:
    
                                        name="%s.%s" % (x,'pdb')
                                        namezipped="%s.%s" % (x,'pdb.zip')
                                        zout=os.path.join(loopfolder,x)
                                        ori=os.path.join(loopfolder,x,name)
                                        orizipped=os.path.join(loopfolder,x,namezipped)

                                        folo=os.path.join(loopfolder,x)
                                        #print ori
                                        #print ori
                                        orinpdb=numpdb.NumPdb(ori, {
                                                        "phi_psi": False,
                                                        "sstruc": False,
                                                        "backbone_only": True,
                                                        "protein_only": False,
                                                        "detect_incomplete": False,
                                                        "configuration": False#,
                                                        #"info": False
                                                        })
                                        ccdic={}
                                        result={}
                                        resultnh={}
                                        #resultnh={}
                                        resultlist=[]
                                        seq1=orinpdb.sequence()
                                        if cor==1:
                                            loopcorrel_name='loop_correl'
                                        else:
                                            loopcorrel_name="%s%s " % ('loop_correl',cor)
                                           # print loopcorrel_name
                                        loops=os.path.join(loopfolder,x,'loop_correl/crosscorrelation/loops_repair')
                                       # print loops
                                        
                                        ccsort=os.path.join(loopfolder,x,loopcorrel_name,'crosscorrelation/ccsort.cpv')
                                        #print cor, ccsort
                                        linkertxt=os.path.join(loopfolder,x,'link_it/edited_linker.txt')
                                        oricc=os.path.join(loopfolder,x,loopcorrel_name,'oricrosscorrelation.cpv')
                                        if os.path.isfile(ccsort):
                                            #print cor
                                            Dssp(ori,output_dir=folo)
                                            oridssp='%s%s' % (ori[:-3],'dssp')
                                            orirec=parse_dssp(oridssp)
                                            with open (ccsort,'r') as hurz:
                                                lines=hurz.readlines()
                                               # print ori
                                                nr=lines[1].split()[0]
                                                corr=lines[1].split()[2]
                                                meancorrel.append(float(corr))
                                                for j in lines[1:]:
                                                    u=int(j.split()[0])
                                                    
                                                    cci=j.split()[2]
                                                    
                                                    ccdic[u]=float(cci)
                                               
                                            with open (linkertxt, 'r') as linki:
                                                lineslinker=linki.readlines()
                                            oriccval=0
                                            h=os.listdir(loops)
                                            #for i in h:
                                            #    if i.endswith ('.zip'):
                                            #        pdbzipname=os.path.join(loops,i)
                                            #        #print pdbzipname
                                            #        with zipfile.ZipFile(pdbzipname,'r') as zi:
                                            #            zi.extractall(loopfolder)
                                            #        os.remove(pdbzipname)
                                            h=os.listdir(loops)
                                            res1=int(x.split('_')[0])
                                            res2=int(x.split('_')[1])
                                           # print res1,res2
                                            for i in h:
                                                if i.endswith('.pdb'):
                                                    #print cor, i
                                                    lino=int(i[0:3])
                                                    lpath=os.path.join(loops,i)
                                                    npdbl = numpdb.NumPdb(lpath,{
                                                                "phi_psi": False,
                                                                "sstruc": False,
                                                                "backbone_only": True,
                                                                "protein_only": False,
                                                                "detect_incomplete": False,
                                                                "configuration": False,
                                                                "info": False
                                                                })
                                                    npdb=npdbl.copy(resno=[res1,res2])
                                                    lnr=int(i[0:3]) *4+1
                                                    seqnr=int(i[0:3]) *4
    
                                                    rootm2=rmsd(orinpdb['xyz'],npdb['xyz'])
    
                                                    try:
    
                                                        tms=tmscore(orinpdb['xyz'],npdb['xyz'])
                                                    except:
                                                        tms=0
                                                        #print l,x,i
                                                      
                                                    seq2=lineslinker[seqnr].strip()
                                                    
                                                    #print seq1
                                                    #sequence identity
                                                    si=0
                                                    for sf in range(0, len(seq1),1):
                                                        if seq1[sf]==seq2[sf+1]:
                                                            si+=1
                                                    sqi=(si/len(seq1))*100
                                                    #print seq1,"and",seq2[1:-1]
                                                    #print sqi
                                                        
                                                    
                                                        
                                                    #print rootm2
                                                    result[i]=[rootm2]
                                                    result[i].append(lineslinker[lnr].strip())
                                                    result[i].append(sqi)
                                                    result[i].append(tms)
                                                    result[i].append(ccdic[lino])
                                                    #result[i].append(oriccval)
                                                    resultlist.append(rootm2)
                                                    
                                                else:
                                                    continue
                                                
                                            #den result dictonary auf eintraege resuzieren die nicht auf der Liste stehen:
                                            resultnh=result.copy()
                                            for key in result:
                                                if result[key][1].upper() in homo:
                                                    del resultnh[key]
                                            #einen dictonary erstellen wo die homologen linker nicht drin sind:
                                            resultnhli=result.copy()
                                            for key in result:
                                                if result[key][2] >=75:
                                                    del resultnhli[key]
                                    
                                            pdbkey="%s_%s.%s" % (nr.zfill(3),'bb','pdb')
    
                    
                                            best=sorted(resultlist)
    
                                            nr_opt=min(result, key=result.get)
                                            
                                            i=1
                                            
                                            # an welcher Stelle is loop nr 1 gerankt
                                            
                                            for c,line in enumerate (lines, 1):
                                                gu=line.split()
                                                if gu[0]==str(nr_opt)[1:3]:
                                                    z=c
                                                    rank.append(c-1)
    
                                                i+=1
    
                                            if pdbkey==min(result, key=result.get):
                                                opti=1
                                            else:
                                                opti=0
                                            if result[pdbkey][0]<result['001_bb.pdb'][0]:
                                                better=1
                                            else:
                                                better=0
                                            if result[pdbkey][0]-result['001_bb.pdb'][0]>4:
                                                cases.write(ori + os.linesep)
                                                #print ori
                                            #print 'min', min(result, key=result.get)
                                            #print 'result',result['001_bb.pdb'][0]
                                            #print 'pdbkey' ,pdbkey
                                            rline="%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % (x,result['001_bb.pdb'][0],pdbkey,result[pdbkey][0],min(result, key=result.get),best[0],z-1,best[-1],opti,better,'\n')
                                            bestrmsds.append(result[pdbkey][0])
                                            seqidentfull.append(result[pdbkey][2])
                                            wodens.append(result['001_bb.pdb'][0])
                                            optimalrmsd.append(best[0])
                                            ###nur extragute
                                            if corrband2<float(corr) <corrband1 and result[pdbkey][2]<75:
                                                extragut.append(result[pdbkey][0])
                                            #
                                            ana.write(rline)
                                            
                                            
                                            ##############
                                            
                                            #schleife ueber die loops
                                            
                                            for y in range (1,100,1):
                                                loop="%s_%s.%s" % (str(y).zfill(3),'bb','pdb')
                                               # print loop
                                                if loop in resultnh:
                                                    wodensnh.append(result[loop][0])
                                                    break
                                                else:
                                                    continue
                                                
                                            # das ganze fuer mit density
                                            
                                            for line in lines[1:]:
                                                #lines.next()
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                
                                               # print 'richtig so?', result[lonam][1].upper()
                                                if result[lonam][1].upper() not in homo:
                                                    #print 'bloo'
                                                    bestrmsdsnh.append(result[lonam][0])
                                                    break
                                                else:
                                                    continue
                                            #homologe linker aussortieren:
                                            
                                            for y in range (1,100,1):
                                                loop="%s_%s.%s" % (str(y).zfill(3),'bb','pdb')
                                                #print loop
                                               # print loop
                                                if loop in resultnhli:
                                                    bestlinkitnohomo=loop
                                                    wodensnhli.append(result[loop][0])
                                                    #print loop
                                                    hfj=os.path.join(loops,loop)
                                                    hjd='%s%s' % (hfj[:-3],'dssp')
                                                    Dssp(hfj,output_dir=loops)
                                                    dsrec=parse_dssp(hjd)
                                                    
                                                  #  print hjd
                                                    
                                                    #with open ('hjd') as ds:
                                                    #    gd=ds.readlines()
                                                    #    for line in gd:
                                                    #        if line.split() []
                                                    if any('H'==e[3]for e in dsrec[1:])==True:
                                                        #print 'found helical loop'
                                                        linkitheli.append(result[loop][0])
                                                        #fragfitheli.append
                                                    else:
                                                        linkitother.append(result[loop][0])
                                                    #for i in dsrec[1:]:
                                                    #    if i[3]=='H':
                                                    #        print i[3]
                                                    #print "dssp", dsrec[1]
                                                    break
                                                else:
                                                    continue
                                            for line in lines[1:]:
                                                #lines.next()
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                hfj2=os.path.join(loops,lonam)
                                                Dssp(hfj2,output_dir=loops)
                                                hjd2='%s%s' % (hfj2[:-3],'dssp')
                                                dsrec2=parse_dssp(hjd2)
                                               # print dsrec2
                                                ccsortbest=result[lonam][0]
                                                if result[lonam][2]<75:                                               
                                                    if any('H'==e[3]for e in dsrec2[1:])==True:
                                                        fragfitheli.append(result[lonam][0])
                                                    elif any ('E'==e[3]for e in dsrec2[1:])==True:
                                                        fragfitbeta.append(result[lonam][0])
                                                    else:
                                                        fragfitother.append(result[lonam][0])
                                                    break
                                                else:
                                                    continue
                                            optrmsdnhnr=min(resultnhli,key=resultnhli.get)
                                            optrmsdnhl=result[optrmsdnhnr][0]
                                           # print optrmsdnhl
                                            for line in lines[1:]:
                                                #lines.next()
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                if result[lonam][2]<75:
                                                    
                                                    bestrmsdnhli.append(result[lonam][0])
                                                    bestfragnoh=result[lonam][0]
                                                    tmsnhli.append(result[lonam][3])
                                                    meancc75.append(result[lonam][4])
                                                    diffh75= abs(float(ccdic[int(no)])-float(oriccval))
                                                    diffh75full=float(ccdic[int(no)])-float(oriccval)
                                                    ordiff75.append( diffh75 )
                                                    ordiffh75full.append(diffh75full )
                                                    suflagnh=0
                                                    if bestlinkitnohomo==optrmsdnhl:
                                                        continue
                                                    else:
                                                        linkitonlybetter.append(result[bestlinkitnohomo][0])
                                                        fragfitonlybetter.append(result[lonam][0])
                                                    if result[lonam][0]==optrmsdnhl:
                                                        sucounternh+=1
                                                        suflagnh=1
                                                    else:
                                                        suflagnh=0   
                                            
                                                
                                                    
                                                    try:
                                                        zhf=0
                                                        if any ('H'==e[3] for e in orirec[1:])==True:
                                                            fragfithelior.append(result[lonam][0])
                                                            zhf=1
                                                        elif any ('E'==e[3] for e in orirec[1:])==True:
                                                            fragfitbetaor.append(result[lonam][0])
                                                            zhr=2
                                                        else:
                                                            fragitotheror.append(result[lonam][0])
                                                            zhf=3
                                                        diffile.write(str(float(ccdic[int(no)])-float(oriccval)) +',')
                                                        sequenceidentity.append(result[lonam][2])
                                                        
                                                        #if math.isnan(diffh75):
                                                        #    print loopfolder
                                                        #    break
                                                        #print diffh75
                                                        break
                                                    except:
                                                        continue
                                                else:
                                                    continue
                                            suflag=0
                                            if result[pdbkey][0] ==best[0]:
                                                suflag=1
                                                sucounter +=1
                                            else:
                                                suflag=0
                                            linkitfirst.append(result['001_bb.pdb'][0])
                                            fragfitfirst.append(ccsortbest)
                                            if result['001_bb.pdb'][0]> 5:
                                                linkitfirstbad.append(result['001_bb.pdb'][0])
                                                fragfitfirstbad.append(ccsortbest)
                                            if result['001_bb.pdb'][0]>rmsdval[str(l)]:
                                                linkitfirstbadf.append(result['001_bb.pdb'][0])
                                                fragfitfirstbadf.append(ccsortbest)
                                                badf+=1
                                            #looline= "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % (g,x,lonam,result[lonam][1],result[lonam][2],result[lonam][0],zhf,str(float(ccdic[int(no)])-float(oriccval)),result[lonam][3],diffh75full,nr_opt,best[0],pdbkey,result[pdbkey][0],suflag,optrmsdnhl,suflagnh,result['001_bb.pdb'][0],ccsortbest,'\n')
                                            #loofile.write(looline)
                                            for line in lines[1:]:
                                                #lines.next()
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                
                                                if result[lonam][2]<90:
                                                    fragfithomo90.append(result[lonam][0])
                                                    #tmsnhli.append(result[lonam][3])
                                                    break
                                                else:
                                                    continue
                                            
                                            for line in lines[1:]:
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                if  int(no)<50 and result[lonam][2]<75:
                                                    result50.append(result[lonam][0])
                                                    break
                                                else:
                                                    continue
                                            
                                            for line in lines[1:]:
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                if  int(no)<20 and result[lonam][2]<75:
                                                    result20.append(result[lonam][0])
                                                    break
                                                else:
                                                    continue
                                            #die besten five alternative
                                            bestfive={}
                                            for line in lines[1:]:
                                            #    #lines.next()
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                if result[lonam][2]<75:
                                                    
                                                    bestfive[lonam]=(result[lonam])
                                                   
                                                    if len(bestfive)==5:
                                                        #print bestfive
                                                        
                                                       # print sorted(bestfive.values())[0] 
                                                        bestfiveall.append(sorted(bestfive.values())[0] [0])
                                                        diffbestfive.append(sorted(bestfive.values())[0] [4])
                                                        break
                                                else:
                                                    continue
                                            #die besten zehn
                                            #bestfive=[]
                                            #for line in lines[1:]:
                                            #    #lines.next()
                                            #    no=line.split()[0]
                                            #    lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                            #    if result[lonam][2]<75:
                                            #        bestfive.append(float(result[lonam][0]))
                                            #        if len(bestfive)==5:
                                            #            #print sorted(bestfive)
                                            #            bestfiveall.append(sorted(bestfive)[0])
                                            #            break
                                            #    else:
                                            #        continue
                                            #die besten zehn
                                            bestten=[]
                                            for line in lines[1:]:
                                                #lines.next()
                                                no=line.split()[0]
                                                lonam="%s_%s.%s" % (no.zfill(3),'bb','pdb')
                                                if result[lonam][2]<75:
                                                    bestten.append(float(result[lonam][0]))
                                                    if len(bestten)==10:
                                                        #print sorted(bestfive)
                                                        besttenall.append(sorted(bestten)[0])
                                                        break
                                                else:
                                                    continue       # print '.'
                                            #helix_sheet_by header
                                            for key in helix:
                                                        if helix[key] [0] >=l and helix[key] [2]==g and helix[key] [3] in range (y1-2,y1+2):
                                                            #print key, helix[key], g, y1
                                                            helixlinkit.append(result['001_bb.pdb'][0])
                                                            helixfragfit.append(ccsortbest)
                                                            
                                            for key in sheet:
                                                if sheet[key] [0] >=l and sheet[key] [2]==g and sheet[key] [3] in range (y1-2,y1+2):
                                                    #print key, sheet[key], g, y1
                                                    sheetlinkit.append(result['001_bb.pdb'][0])
                                                    sheetfragfit.append(ccsortbest)
                                            for key in fullsheet:
                                                if fullsheet[key] [3] >=l and fullsheet[key] [0]==g and fullsheet[key] [1] in range (y1-2,y1+2):
                                                    #print key, sheet[key], g, y1
                                                    fullsheetlinkit.append(result['001_bb.pdb'][0])
                                                    fullsheetfragfit.append(ccsortbest)
                    #print 'hallo'
                    if len(linkitheli)==0:
                        linkitheli.append(0)
                    if len(extragut)==0:
                        extragut.append(0)
                    if len(bestrmsds)==0:
                        bestrmsds.append(0)
                    if  len(wodens)==0:
                        wodens.append(0)
                    if len(optimalrmsd)==0:
                        optimalrmsd.append(0)
                    if len(bestrmsdsnh)==0:
                        bestrmsdsnh.append(0)
                    if len(wodensnh)==0:
                        wodensnh.append(0)
                    if len(bestrmsdnhli)==0:
                        bestrmsdnhli.append(0)
                    if len(wodensnhli)==0:
                        wodensnhli.append(0)
                    if  len(meancorrel)==0:
                        meancorrel.append(0)
                    if len(bestfiveall)==0:
                        bestfiveall.append(0)
                    if len(linkitother)==0:
                        linkitother.append(0)
                    if len(rank)==0:
                        rank.append(0)
                    if len(tmsnhli)==0:
                        tmsnhli.append(0)
                    if len(result50)==0:
                        result50.append(0)
                    if len(fragfithomo90)==0:
                        fragfithomo90.append(0)
                    if len (result20)==0:
                        result20.append(0)
                    if len (besttenall)==0:
                        besttenall.append(0)
                    if len (ori_meancorrel)==0:
                        ori_meancorrel.append(0)
                    if len (or_diffcc)==0:
                        or_diffcc.append(0)
                    if len (ordiff75)==0:
                       # print loopfolder
                        ordiff75.append(0)
                    if len ( ordiffh75full)==0:
                         ordiffh75full.append(0)
                    if len (diffbestfive)==0:
                        diffbestfive.append(0)
                    if len (meancc75)==0:
                        meancc75.append(0)
                    if len(fragfitheli)==0:
                        fragfitheli.append(0)
                    if len(fragfitother)==0:
                        fragfitother.append(0)
                    if len(sequenceidentity)==0:
                        sequenceidentity.append(0)
                    if len(seqidentfull)==0:
                        seqidentfull.append(0)
                    if len(fragfitbeta)==0:
                        fragfitbeta.append(0)
                    if len(fragfithelior)==0:
                        fragfithelior.append(0)
                    if len(fragfitbetaor)==0:
                        fragfitbetaor.append(0)
                    if len(fragitotheror)==0:
                        fragitotheror.append(0)
                    if len(linkitonlybetter)==0:
                        linkitonlybetter.append(0)
                    if len(fragfitonlybetter)==0:
                        fragfitonlybetter.append(0)
                    if len(fragfitfirst)==0:
                        fragfitfirst.append(0)
                    if len(linkitfirst)==0:
                        linkitfirst.append(0)
                    if len(linkitfirstbad)==0:
                        linkitfirstbad.append(0)    
                    if  len(fragfitfirstbad)==0:
                        fragfitfirstbad.append(0)
                    if len(linkitfirstbadf)==0:
                        linkitfirstbadf.append(0)
                    if len(fragfitfirstbadf)==0:
                        fragfitfirstbadf.append(0)
                        
                    if len(helixlinkit)==0:
                        helixlinkit.append(0)
                    if len(helixfragfit)==0:
                        helixfragfit.append(0)
                    if len(sheetlinkit)==0:
                        sheetlinkit.append(0)
                    if len(sheetfragfit)==0:
                        sheetfragfit.append(0)
                    if len(fullsheetlinkit)==0:
                        fullsheetlinkit.append(0)
                    if len(fullsheetfragfit)==0:
                        fullsheetfragfit.append(0)
                    st,sp = scipy.stats.ttest_rel(wodensnh,bestrmsdsnh)
                    print st,sp
                    break
                    #print fragfitonlybetter,linkitonlybetter
                    #break
                   # print 'looplaenge:' ,l
                   # print 'summe ordiff75:', sum(ordiff75)
                   # print 'length ordiff75:', len(ordiff75)
                   # print 'geteilt:',sum(ordiff75)/len(ordiff75)
                    loopline= "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s" % (l,sum(bestrmsds)/len(bestrmsds),np.std(bestrmsds),
                                                             sum(wodens)/len(wodens),np.std(wodens),
                                                             sum(bestrmsdsnh)/len(bestrmsdsnh),np.std(bestrmsdsnh),
                                                             sum(wodensnh)/len(wodensnh),np.std(wodensnh),
                                                             sum(bestrmsdnhli)/len(bestrmsdnhli),sum(fragfithomo90)/len(fragfithomo90),np.std(bestrmsdnhli),
                                                             sum(wodensnhli)/len(wodensnhli),np.std(wodensnhli),
                                                             len(wodensnh),sum (optimalrmsd)/len(optimalrmsd),sum(meancorrel)/len(meancorrel),
                                                             sum(extragut)/len(extragut),sum(bestfiveall)/len(bestfiveall),
                                                             sum(besttenall)/len(besttenall),sum(linkitheli)/len(linkitheli),
                                                             sum(linkitother)/len(linkitother),sum(rank)/len(rank),
                                                             sorted(bestrmsdnhli)[len(bestrmsdnhli)//2],sum(tmsnhli)/len(tmsnhli),
                                                             sum(result50)/len(result50),sum(result20)/len(result20),sum(ori_meancorrel)/len(ori_meancorrel),
                                                             sum(or_diffcc)/len(or_diffcc),sum(ordiff75)/len(ordiff75),sum(ordiffh75full)/len(ordiffh75full),sum(diffbestfive)/len(diffbestfive),
                                                             sum(meancc75)/len(meancc75),sum(fragfitheli)/len(fragfitheli),sum(fragfitother)/len(fragfitother),sum(fragfitbeta)/len(fragfitbeta),
                                                             sum(sequenceidentity)/len(sequenceidentity),sum(seqidentfull)/len(seqidentfull),sum(fragfithelior)/len(fragfithelior),sum(fragfitbetaor)/len(fragfitbetaor),
                                                             sum(fragitotheror)/len(fragitotheror),sucounter/len(bestrmsds),sucounternh/len(bestrmsds),sum(linkitonlybetter)/len(linkitonlybetter),
                                                             sum(fragfitonlybetter)/len(fragfitonlybetter),sum(linkitfirst)/len(linkitfirst),sum(fragfitfirst)/len(fragfitfirst),
                                                             sum (linkitfirstbad)/len(linkitfirstbad),sum(fragfitfirstbad)/len(fragfitfirstbad),sum (linkitfirstbadf)/len(linkitfirstbadf),sum(fragfitfirstbadf)/len(fragfitfirstbadf),badf,
                                                             sum(helixlinkit)/len(helixlinkit),sum(helixfragfit)/len(helixfragfit),len(helixfragfit),sum(sheetlinkit)/len(sheetlinkit),sum(sheetfragfit)/len(sheetfragfit),len(sheetfragfit),
                                                             sum(fullsheetlinkit)/len(fullsheetlinkit),sum(fullsheetfragfit)/len(fullsheetfragfit),len(fullsheetfragfit),'\n')
                    #print np.std(bestrmsdsnh) 
                    analyse.write(loopline)
                    analyse.close
                    #linkertxt.close()
                    diffile.write('\n')
                    valu.write(str(l))
                    valu.write(',')
                    valu.write('\n')
                    for j in (bestrmsdnhli):
                        valu.write(str(j))
                        valu.write(',')
                        valu.write('\n')
