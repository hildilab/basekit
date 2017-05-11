import sys
import auto_parameterization as ap

def main(input_file, outdir):
    ligands,aligands=ap.get_ligands(input_file ,outdir)
    if len(ligands)>0:
        ap.main(input_file ,outdir )
        with open(outdir+ '/Ligands.txt' ,'r') as f:
            with open(outdir+ '/files/atomdict.db' ,'a') as o:
                for line in f:
                    with open(outdir+'/'+ line.split()[0]+'.db' ,'r') as ff:
                        for line2 in ff:
                            o.write(line2)
        input_file="Dowser++_"+input_file
    print(input_file)

if __name__ == "__main__":
    if len(sys.argv)!=3:
        print("Error! Too many or too less arguments! Two files are needed.")
    else:
        main(sys.argv[1],sys.argv[2])
