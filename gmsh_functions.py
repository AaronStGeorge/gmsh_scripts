import dolfin as dolf
import pylab as p
import subprocess as s
import sys, os

def pos_from_xml(mesh,results):
    
    # read files into dolfin format
    mesh = dolf.Mesh(mesh)
    Q    = dolf.FunctionSpace(mesh,'CG',1)
    field= dolf.Function(Q)
    dolf.File(results) >> field
    results = results.rstrip("xml")+"pos"
    output = open(results,'w')
    cell_type = mesh.type().cell_type()
    nodes = mesh.coordinates()
    n_nodes = mesh.num_vertices()
    nodes = p.hstack((nodes,p.zeros((n_nodes,3 - p.shape(mesh.coordinates())[1]))))
    f = lambda a: [field(a[0],a[1])]
    nodes = p.hstack((nodes,map(f,nodes)))
    cells = mesh.cells()
    n_cells = mesh.num_cells()

    # write each triangular elment and the values at each corner   
    output.write("View \"NodalValues\" { \n")
    for ii,cell in enumerate(cells):

        # ST(first vierice of triangual elment, second, third){value at first, second, third};\n
        output.write("ST({0:g},{1:g},{2:g},{3:g},{4:g},{5:g},{6:g},{7:g},{8:g}){{{9:g},{10:g},{11:g}}};\n".format(nodes[cell[0]][0],  nodes[cell[0]][1],  nodes[cell[0]][2],  nodes[cell[1]][0],  nodes[cell[1]][1],  nodes[cell[1]][2],  nodes[cell[2]][0],  nodes[cell[2]][1],  nodes[cell[2]][2],  nodes[cell[0]][3],  nodes[cell[1]][3],  nodes[cell[2]][3]))
   
    output.write("};")

    output.close()
    
    # to get the rounding on the numbers write we mus oepn it and re-save in gmsh
    output = open('temp.geo','w')
    output.write("Include '{0}'; \n".format(results))
    output.close()
    print "please click on post processing view and hit save"
    s.call(['gmsh','temp.geo'])
    os.remove('temp.geo')

def xml_to_gmsh(mesh,results=None):
    
    """
    function iterates through mesh and results (field) data both provided in .xml format and writes a output .msh file with the field data read in as node data. The output file same name as the field file passed or if no filed file passed will have name of input mesh.

    :param mesh:  (.xml) Mesh that is to be written to a file
    :param field: (.xml) results data
    """
    
    # get output file
    if results == None:
        fname = mesh.rstrip("xml")+"msh" 
        output = open(fname,'w')    
    else:
        fname = results.rstrip("xml")+"msh" 
        output = open(fname,'w')
    
    # read files into dolfin format
    mesh = dolf.Mesh(mesh)
    Q    = dolf.FunctionSpace(mesh,'CG',1)
    if results != None:
        field= dolf.Function(Q)
        dolf.File(results) >> field

    cell_type = mesh.type().cell_type()

    nodes = mesh.coordinates()
    n_nodes = mesh.num_vertices()

    nodes = p.hstack((nodes,p.zeros((n_nodes,3 - p.shape(mesh.coordinates())[1]))))

    cells = mesh.cells()
    n_cells = mesh.num_cells()

    output.write("$MeshFormat\n" + 
                "2.2 0 8\n" +
                "$EndMeshFormat\n" +
                "$Nodes \n" +
                "{0:d}\n".format(n_nodes))

    for ii,node in enumerate(nodes):
        output.write("{0:d} {1} {2} {3}\n".format(ii+1,node[0],node[1],node[2]))

    output.write("$EndNodes\n")

    output.write("$Elements\n" + 
    "{0:d}\n".format(n_cells))
    
    for ii,cell in enumerate(cells):

        if cell_type == 1:
            output.write("{0:d} 1 0 {1:d} {2:d}\n".format(ii+1,int(cell[0]+1),int(cell[1]+1)))

        elif cell_type == 2:
            output.write("{0:d} 2 0 {1:d} {2:d} {3:d}\n".format(ii+1,int(cell[0]+1),int(cell[1]+1),int(cell[2]+1)))

        elif cell_type == 3:
            output.write("{0:d} 4 0 {1:d} {2:d} {3:d} {4:d}\n".format(ii+1,int(cell[0]+1),int(cell[1]+1),int(cell[2]+1),int(cell[3]+1)))

        else:
            print "Unknown cell type"

    output.write("$EndElements\n")

    if results != None: 
        output.write("$NodeData\n" +
                    "1\n" +
                    "\"NodalValues\"\n" +
                    "0\n" + # zero real tags 
                    "3\n" + # three integer tags
                    "0\n" + # the time step
                    "1\n" + #1-component (scalar) field
                    "{0:d}\n".format(len(nodes)))
        for ii,node in enumerate(nodes):
            output.write("{0:d} {1:g}\n".format(ii+1,field(node[0],node[1])))

        output.write("$EndNodeData\n")    

    output.close()
    
    # There is some numerical precision error that prevents files created by
    # this script from being converted back into xml by dolfin-convert opening the file in gmsh and resaving fixes this 
    print 'calling gmsh...'
    s.call(['gmsh',fname,'-2'])
    
def hess_pos(mesh,results):

    # check if results is .xml or .pos
    if results[-3:] == 'xml':
        print 'results in .xml format searching for .pos format...'
        try:
            with open(results[:-3]+'pos'):
                print 'found...'
                pass
        except IOError:
            print 'none found creating .pos file...'
            if mesh[-3:] == 'msh':
                print 'searching for mesh file in .xml format in cwd (needed to convert results to .pos)...'
                try:
                    with open(mesh[:-3]+'xml'): 
                        'found...  converting'
                        pos_from_xml(mesh[:-3]+'xml',results)
                except IOError:
                    'none found creating .xml file...'
                    s.call(["dolfin-convert", mesh, mesh[:-3]+'xml'])   
                    'converting...'
            pos_from_xml(mesh[:-3]+'xml',results)
        results = results[:-3] + 'pos'          

    # check if mesh is .xml or .msh
    if mesh[-3:] == 'xml':
        print 'mesh in .xml format searching for .msh format...'
        try:
            with open(mesh[:-3]+'msh'):
                print 'found...'
        except IOError:
                print 'none found, converting...'
                xml_to_gmsh(mesh)
        mesh = mesh[:-3] + 'msh'
    
    # write .geo file to create .pos files
    output = open('create_hess_pos.geo','w')
    output.write( "Merge '{0}'; \n".format(mesh))
    output.write( "Include '{0}'; \n".format(results))
    output.write( "Background Mesh View[0]; \n")
    output.write( "Field[1].CropNegativeValues = 0; \n")
    output.write( "Field[2] = MaxEigenHessian; \n")
    output.write( "Field[2].Delta = 0.0001; \n")
    output.write( "Field[3] = MathEval; \n")
    output.write( "Field[3].F = \"abs(F2) + 1\"; \n")
    output.write( "Field[4] = MathEval; \n")
    output.write( "Field[4].F = \"log(F3) \"; \n")
    output.close()
    
    # run .geo file
    # this badly needs to be removed... there must be some way create and save a
    #   view from a .geo file
    print "sorry to do this... please hit \"m\", then in the mesh pane hit define>size_fields. A window should pop choose from the list the last MathEval field. Then in the lower left hit view then create new view. This will take a few minutes. Once it is done in the post-processing pane there will be a new view click on this and select \"save as\" copy paste \"Thresh_on_log_abs_MaxEigenHessian_plus_one.pos \" (without the quotes) into the save as dialog box, then hit save."
    s.call(["gmsh","create_hess_pos.geo"])



