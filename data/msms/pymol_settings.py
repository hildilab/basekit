import pymol
import sys
sys.path.append( "${TMPL_DIR}" )
import wfmesh
pymol.finish_launching()
pymol.cmd.bg_color(color="white")
neighbours=${neighbours}
mean_dic=${mean_dct}
mean_list=${mean_list}
obj_list=${obj_list}
pymol.cmd.load("${nh_file}")
name="${nh_file}".split('/')[-1][:-4]
print name, "${nh_file}"
if "${pref}"!="":
    newname="${pref}_"+name
    pymol.cmd.set_name(name, newname)
pymol.cmd.hide( representation="line")
pymol.cmd.show( representation="cartoon")

for index, elem in enumerate(mean_dic):
    pymol.cmd.select( "temp", "(id "+str(elem)+" and resname NEH)" )
    pymol.cmd.alter("temp", "vdw="+str(mean_dic[elem]))
    pymol.cmd.rebuild()
for index, elem in enumerate(neighbours):
    liste=""
    for neighbour in neighbours[elem]:
        if liste!="":
            liste=liste+"+"
        liste=liste+"(id "+str(neighbour[0])+" and resi "+str(neighbour[1])+" and chain "+neighbour[2]+")"
    if "${pref}"!="":
        pymol.cmd.select( "${pref}_nb_"+str(mean_list[index]), liste )
    else:
        pymol.cmd.select( "neighbour_"+str(mean_list[index]), liste )
for index, elem in enumerate(mean_list):
    atomno, resno, chain = elem.split('_')
    liste="(resi "+str(resno)+" and resname NEH)"
    if "${pref}"!="":
        pymol.cmd.create("${pref}_neh_"+str(mean_list[index]), liste)
        pymol.cmd.show(representation="spheres", selection="${pref}_neh_"+str(mean_list[index]))
        pymol.cmd.set("sphere_transparency", value=0.5, selection="${pref}_neh_"+str(mean_list[index]))
    else:
        pymol.cmd.create("neh_"+str(mean_list[index]), liste)
        pymol.cmd.show(representation="spheres", selection="neh_"+str(mean_list[index]))
        pymol.cmd.set("sphere_transparency", value=0.5, selection="neh_"+str(mean_list[index]))
if "${pref}"!="":
    pymol.cmd.select( "${pref}_neh", "${pref}_neh_*")
else:
    pymol.cmd.select( "neh", "neh_*")
for index, elem in obj_list:
    if "${pref}"!="":
        wfmesh.createWFObj(elem,"${pref}_obj_"+str(mean_list[index-1]))
    else:
        wfmesh.createWFObj(elem,"obj_"+str(mean_list[index-1]))
pymol.cmd.delete("temp")
if "${pref}"!="":
    pymol.cmd.group("${pref}_nb", "${pref}_nb_*")
    pymol.cmd.group("${pref}_obj", "${pref}_obj_*")
    pymol.cmd.group("${pref}_nehs", "${pref}_neh_*")
else:
    pymol.cmd.group("neighbours", "neighbour_*")
    pymol.cmd.group("obj", "obj_*")
    pymol.cmd.group("neighbourholes", "neh_*")
pymol.cmd.select("resname NEH")
pymol.cmd.order("*", "yes")
if "${pref}"!="":
    pymol.cmd.order("order *_"+name+" sele *_neh *_nehs *_nb *_obj")
    pymol.cmd.center(newname)
    pymol.cmd.zoom(newname)
else:
    pymol.cmd.order("order ${nh_file} sele neh neighbourholes neighbours obj")
    pymol.cmd.center(name)
    pymol.cmd.zoom(name)
pymol.cmd.spectrum( "b", "blue_white_red", "sele")
