import pymol
import sys
sys.path.append( "${TMPL_DIR}" )
import wfmesh
pymol.finish_launching()
pymol.cmd.bg_color(color="white")
neighbours=${neighbours}
mean_dic=${mean_dct}
obj_list=${obj_list}
pymol.cmd.load("${nh_file}")
pymol.cmd.hide( representation="line")
pymol.cmd.show( representation="cartoon")
pymol.cmd.select( "neh", "resname NEH")
pymol.cmd.show(representation="spheres", selection="neh")
pymol.cmd.spectrum( "b", "blue_white_red", "neh")
for index, elem in enumerate(mean_dic):
    pymol.cmd.select( "temp", "(resi "+str(index)+" and resname NEH)" )
    pymol.cmd.alter("temp", "vdw="+str(mean_dic[elem]))
    pymol.cmd.rebuild()
for index, elem in enumerate(neighbours):
    liste=""
    for neighbour in neighbours[elem]:
        if liste!="":
            liste=liste+"+"
        liste=liste+"(id "+str(neighbour[0])+" and resi "+str(neighbour[1])+" and chain "+neighbour[2]+")"
    pymol.cmd.select( "neighbour_"+str(elem), liste )
for index, elem in obj_list:
    wfmesh.createWFObj(elem,"obj_"+str(index))
pymol.cmd.delete("temp")
pymol.cmd.group("neighbours", "neighbour_*")
pymol.cmd.group("obj", "obj_*")
pymol.cmd.select("resname NEH")
pymol.cmd.order("*", "yes")
pymol.cmd.order("order ${nh_file} sele neh neighbours obj")
pymol.cmd.center("${nh_file}".split('/')[-1][:-4])
pymol.cmd.zoom("${nh_file}".split('/')[-1][:-4])
