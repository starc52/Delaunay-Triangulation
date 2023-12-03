from dcel.dcel import Dcel
from dcel.dcel import *

def edge_swap(dcel):
    num_hedges = len(dcel.hedges)
    iteration = 0
    iterator = -1
    # iterate until convergence
    while True:
        iteration += 1
        iterator += 1
        iterator %= num_hedges
        # iterate over half-edges
        hedge = dcel.hedges[iterator]
        # find a non-legal edge
        if hedge.external or hedge.legal:
            continue
        else:
            # implement the incircle test
            v3 = hedge.origin
            v1 = hedge.nexthedge.origin
            v2 = hedge.prevhedge.origin
            q = hedge.twin.prevhedge.origin
            # test legality
            if dcel.incircle(v1, v2, v3, q) == 1:
                dcel.delete_edge(v1, v3)
                dcel.add_edge(v2, q)
            else:
                hedge.legal = True
        # self.saveplot(f"ex1_{iteration}.eps")
        # break testing
        count_hedges = 0
        for hedge in dcel.hedges:
            if hedge.legal == True:
                count_hedges += 1
        if count_hedges == num_hedges:
            break
    print(f"Edge Swap Algorithm Finished in {iteration} iterations over")


def test():
    d = Dcel()
    d.load('ex1.txt')
    d.saveplot('ex1.eps')
    faces_wo_external = [face for face in d.faces if not face.external]
    for a, p, face in zip(d.areas(), d.perimeters(), faces_wo_external):
        print(a, p, [(v.x, v.y) for v in face.vertexlist()])
    d.lawson_delaunay()
    d.saveplot('output.eps')
    faces_wo_external = [face for face in d.faces if not face.external]
    for a, p, face in zip(d.areas(), d.perimeters(), faces_wo_external):
        print(a, p, [(v.x, v.y) for v in face.vertexlist()])


if __name__ == "__main__":
    test()
