using QPSReader, QuadraticModels, QuadraticModelsCPLEX

# path_pb_lp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/netlib"
# path_pb_qp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/marosmeszaros"
path_pb_lp = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib"
# path_pb_qp = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_marosmeszaros"
path_pb_lp_ps = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib_ps"

function presolve_problems_cplex(path_pb, path_pb_ps)
    problems = []
    i_max = 11
    i = 1
    for file_name in readdir(path_pb)
        if file_name[end-3:end] == ".SIF" && !(file_name in ["80BAU3B.SIF" ; "BORE3D.SIF";
                                                        "CAPRI.SIF"; "CZPROB.SIF";
                                                        "ETAMACRO.SIF"; "FINNIS.SIF";
                                                        "FORPLAN.SIF"; "GREENBEA.SIF";
                                                        "GREENBEB.SIF"; "MAROS.SIF";
                                                        "NESM.SIF"; "PEROLD.SIF";
                                                         "PILOT-JA.SIF"; "PILOT-WE.SIF";
                                                         "PILOT.SIF"; "PILOT4.SIF";
                                                         "PILOT87.SIF"; "PILOTNOV.SIF";
                                                          "RECIPELP.SIF"; "SHELL.SIF";
                                                         "SIERRA.SIF"; "STAIR.SIF";
                                                         "STANDATA.SIF"; "STANDGUB.SIF";
                                                        "STANDMPS.SIF"; "TUFF.SIF";
                                                        "VTP-BASE.SIF"; "DTOC3.SIF";
                                                         "HS35MOD.SIF";"QBORE3D.SIF";
                                                        "QCAPRI.SIF"; "QETAMACR.SIF";
                                                          "QFORPLAN.SIF"; "QPCSTAIR.SIF";
                                                        "QPCSTAIR.SIF"; "QPILOTNO.SIF";
                                                        "QRECIPE.SIF"; "QSHELL.SIF";
                                                        "QSIERRA.SIF"; "QSTAIR.SIF";
                                                        "QSTANDAT.SIF"; "UBH1.SIF";
                                                        "YAO.SIF"]) # problems with fixed variables


            println(file_name)
            pb_i = string(path_pb, "/", file_name)
            if file_name in ["BLEND.SIF"; "DFL001.SIF"; "FORPLAN.SIF"; "GFRD-PNC.SIF"; "SIERRA.SIF";
                        "EXDATA.SIF"; "QFORPLAN.SIF"; "QGFRDXPN.SIF"; "VALUES.SIF"]
                qpdata_i = readqps(pb_i, mpsformat=:fixed)
            else
                qpdata_i = readqps(pb_i)
            end
            pb_ps_i = string(path_pb_ps, "/", file_name[1:end-4], "_PS.mps")
            presolve_to_file(QuadraticModel(qpdata_i), path = pb_ps_i)

            if i==i_max
                break
            end
            i += 1
        end
    end
end

presolve_problems_cplex(path_pb_lp, path_pb_lp_ps)
