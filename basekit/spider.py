from __future__ import with_statement
from __future__ import division

import os
import json
import itertools

from utils import copy_dict
from utils.tool import _, _dir_init, PyTool, CmdTool, ScriptMixin
from utils.numpdb import NumPdb
from utils.mrc import get_mrc_header, getMrc

from pdb import *#PdbSplit
#from pdb import PdbEdit, NumpdbTest

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "spider" )
SPIDER_CMD = "spider" 



            
# 2010 Cryo-EM Modeling Challenge: http://ncmi.bcm.edu/challenge


class MrcHeaderPrinter( PyTool ):
    args = [
        _( "mrc_file", type="file", ext="mrc" )
    ]
    def func( self, *args, **kwargs ):
        header = mrc_header( self.mrc_file )
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
        boxsize=getMrc(self.mrc_file,'nx' )
        originx=abs(getMrc(self.mrc_file,'nxstart'))
        originy=abs(getMrc(self.mrc_file,'nystart'))
        originz=abs(getMrc(self.mrc_file,'nzstart'))
        print originx
        size=getMrc(self.mrc_file,'xlen' )
        order=getMrc(self.mrc_file,'mapc' )
        pixelsize=(size/boxsize)
        self._make_script_file(    
            mrc_file=self.relpath( self.mrc_file ),
            order=order,
            pixelsize=pixelsize,
            boxsize=boxsize,
            originx=originx,
            originy=originy,
            originz=originz
        )

        shx = (originx -(boxsize/2)) * pixelsize
        shy = (originy -(boxsize/2)) * pixelsize
        shz = (originz -(boxsize/2)) * pixelsize
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
            bbs=getMrc(self.mrc_file, 'nx')#mrc_dic["nx"]
            obs = bs*ps
            x =  (ol1 - (bbs/2)) * ps - 1 
            y =  (ol2 - (bbs/2)) * ps - 1
            z =  (ol3 - (bbs/2)) * ps - 1
            print [x,y,z, bs,bbs]
        PdbEdit (self.pdb_file, box = [ x, y, z, obs, obs, obs ] )
        

    
class SpiderDeleteFilledDensities( Spider ):
    args = [
        _( "mrc_file", type="file", ext="mrc" ),
        _( "map_file", type="file", ext="cpv" , help= "boxil.cpv" ),
        _( "pdb_file", type="file", ext="pdb" ),
        _( "result_file", type="file" , ext="cpv", help="ergebnisse.cpv"),
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
        boxsize=getMrc(self.mrc_file,'nx' )
        size=getMrc(self.mrc_file,'xlen' )
        pixelsize=(size/boxsize)
        print self.res1
        LoopDelete(self.pdb_file,self.res1['chain'],self.res1['resno']+1,self.res2['resno']-1)
        self._make_script_file( 
            map_name=self.relpath( self.map_file, no_ext=True ), 
            pdb_file=self.relpath( 'noloop.pdb' ),
            result_file=self.relpath( self.result_file, no_ext=True ),
            pixelsize=pixelsize,
            resolution=self.resolution,
            boxsize=boxsize,
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
        self.spider_shift = SpiderShift(
            self.mrc_file,self.pdb_file,
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
        self.pdb_box =  SpiderPdbBox(
           self.spider_shift.edited_pdb_file, 
           self.spider_box.box_file,
           self.spider_shift.map_shift,
           **copy_dict( kwargs, run=False, output_dir=self.subdir("pdbbox"))
        )
        self.spider_delete_filled_densities = SpiderDeleteFilledDensities( 
            self.spider_shift.map_shift,self.spider_box.box_map_file, self.pdb_box.edited_pdb_file, self.spider_box.box_file,
            self.resolution,self.res1,self.res2,
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
        shiftpath=self.relpath(self.spider_shift.output_dir)
        shiftfile="%s/%s" % (shiftpath,'shift.cpv')
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
                        x="%4.3f" % (float(line[31:38])+shxb)
                        a="%7s" % x
                        
                        y="%4.3f" % (float(line[39:46])+shyb)
                        b="%7s" % y
                        z="%4.3f" % (float(line[47:54])+shzb)
                        c="%7s" % z
        
                        line = line = line[0:30] + a + line[38:]
                        line = line = line[0:38] + b + line[46:]
                        line = line = line[0:46] + c + line[53:]
        
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
                                x="%4.3f" % (float(line[31:38])+shxb)
                                a="%7s" % x
                                
                                y="%4.3f" % (float(line[39:46])+shyb)
                                b="%7s" % y
                                z="%4.3f" % (float(line[47:54])+shzb)
                                c="%7s" % z
                        
                                line = line = line[0:30] + a + line[38:]
                                line = line = line[0:38] + b + line[46:]
                                line = line = line[0:46] + c + line[53:]
                        
                                of.write( line )
                        of.write("END")
 