from dcel.dcel import Dcel


def test():
    d = Dcel()
    d.load('ex1.txt')

    if d.hedges[5].external == False and d.hedges[5].twin.external == False:
        v1 = d.hedges[5].origin
        v2 = d.hedges[5].nexthedge.origin
        d.delete_edge(v1, v2)
    d.saveplot('ex1_deleted.eps')


    faces_wo_external = [face for face in d.faces if not face.external]

    for a, p, face in zip(d.areas(), d.perimeters(), faces_wo_external):
        print(a, p, [(v.x, v.y) for v in face.vertexlist()])

    d.add_edge(v1, v2)
    d.saveplot('ex1_added.eps')

    faces_wo_external = [face for face in d.faces if not face.external]

    for a, p, face in zip(d.areas(), d.perimeters(), faces_wo_external):
        print(a, p, [(v.x, v.y) for v in face.vertexlist()])

    d.lawson_delaunay()

if __name__ == "__main__":
    test()
