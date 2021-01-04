using GridapGmsh
using GridapGmsh: gmsh, GmshDiscreteModel

function add_pin!(gmsh, o, r, t, lc)
    factory = gmsh.model.geo

    # inner and outer circle
    l1 = add_circle!(gmsh, o, r, lc)
    l2 = add_circle!(gmsh, o, r + t, lc)

    s1 = factory.addPlaneSurface([l1])
    s2 = factory.addPlaneSurface([l2, l1]) # l1 is a hole

    return s1, s2
end

function add_circle!(gmsh, o, r, lc)
    factory = gmsh.model.geo
    ox, oy = o

    p1 = factory.addPoint(ox, oy, 0, lc)  # center - origin
    p2 = factory.addPoint(ox + r, oy, 0, lc)  # right
    p3 = factory.addPoint(ox, oy + r, 0, lc)  # up
    p4 = factory.addPoint(ox - r, oy, 0, lc)  # left
    p5 = factory.addPoint(ox, oy - r, 0, lc)  # down

    c1 = factory.addCircleArc(p2, p1, p3)
    c2 = factory.addCircleArc(p3, p1, p4)
    c3 = factory.addCircleArc(p4, p1, p5)
    c4 = factory.addCircleArc(p5, p1, p2)

    l1 = factory.addCurveLoop([c1, c2, c3, c4])

    return l1
end

function add_square!(gmsh, s)
    factory = gmsh.model.geo

    p1 = factory.addPoint(0, 0, 0, lc)
    p2 = factory.addPoint(s, 0, 0, lc)
    p3 = factory.addPoint(s, s, 0, lc)
    p4 = factory.addPoint(0, s, 0, lc)

    l1 = factory.addLine(p1, p2)
    l2 = factory.addLine(p2, p3)
    l3 = factory.addLine(p3, p4)
    l4 = factory.addLine(p4, p1)

    cl1 = factory.addCurveLoop([l1, l2, l3, l4])

    return cl1
end

# in cm
p = 1.6        # pitch
p_2 = p / 2    # pitch / 2
ri = 0.5       # internal radius
t = 0.1        # wall thickness
ro = ri + t    # external radius
lc = 0.1

gmsh.initialize()
gmsh.model.add("bwr")
pinTags = Int32[]
cladTags = Int32[]
for i in 1:4, j in 1:4
    if (i, j) == (2, 3) || (i, j) == (3, 2)
        continue
    end
    xo = p_2 + (i - 1) * p
    yo = p_2 + (j - 1) * p
    pinTag, cladTag = add_pin!(gmsh, (xo, yo), ri, t, lc)
    push!(pinTags, pinTag)
    push!(cladTags, cladTag)
end

gdTags = Int32[]
for (i, j) in ((2, 3), (3, 2))
    xo = p_2 + (i - 1) * p
    yo = p_2 + (j - 1) * p
    gdTag, cladTag = add_pin!(gmsh, (xo, yo), ri, t, lc)
    push!(gdTags, gdTag)
    push!(cladTags, cladTag)
end

# h2oTag = add_reflector!(gmsh, 4p)
s = 4p
factory = gmsh.model.geo

p1 = factory.addPoint(0, 0, 0, lc)
p2 = factory.addPoint(s, 0, 0, lc)
p3 = factory.addPoint(s, s, 0, lc)
p4 = factory.addPoint(0, s, 0, lc)

l1 = factory.addLine(p1, p2)
l2 = factory.addLine(p2, p3)
l3 = factory.addLine(p3, p4)
l4 = factory.addLine(p4, p1)

cl1 = factory.addCurveLoop([l1, l2, l3, l4])

h2oTag = factory.addPlaneSurface(vcat(cl1, pinTags, cladTags, gdTags))



pg1 = gmsh.model.geo.addPhysicalGroup(2, pinTags) #! change "pin" to "fuel"
pg2 = gmsh.model.geo.addPhysicalGroup(2, cladTags)
pg3 = gmsh.model.geo.addPhysicalGroup(2, gdTags)
pg4 = gmsh.model.geo.addPhysicalGroup(2, [h2oTag])

gmsh.model.setPhysicalName(2, pg1, "pin") #! change "pin" to "fuel"
gmsh.model.setPhysicalName(2, pg2, "cladding")
gmsh.model.setPhysicalName(2, pg3, "poison")
gmsh.model.setPhysicalName(2, pg4, "water")


# removemos puntos de la geometria duplicados (origenes por ejemplo)
gmsh.model.geo.removeAllDuplicates()

gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(2)

gmsh.write("bwr.msh")






# square cell
factory.addPoint(0, 0, 0, lc, 10)
factory.addPoint(p, 0, 0, lc, 11)
factory.addPoint(p, p, 0, lc, 12)
factory.addPoint(0, p, 0, lc, 13)
factory.addLine(10, 11,  9)
factory.addLine(11, 12, 10)
factory.addLine(12, 13, 11)
factory.addLine(13, 10, 12)
factory.addCurveLoop([9, 10, 11, 12], 3)

factory.addPlaneSurface([1], 1)
factory.addPlaneSurface([2, 1], 2) # le agrego el hole `1` (no estoy teniendo en cuenta nada de sentidos de giro, no se si hay que hacerlo)
factory.addPlaneSurface([3, 2, 1], 3)

# materials
factory.addPhysicalGroup(2, [1], 1)
factory.addPhysicalGroup(2, [2], 2)
factory.addPhysicalGroup(2, [3], 3)
GridapGmsh.gmsh.model.setPhysicalName(2, 1, "pin")
GridapGmsh.gmsh.model.setPhysicalName(2, 2, "cladding")
GridapGmsh.gmsh.model.setPhysicalName(2, 3, "water")

# boundaries
factory.addPhysicalGroup(1,  [9], 4)
factory.addPhysicalGroup(1, [10], 5)
factory.addPhysicalGroup(1, [11], 6)
factory.addPhysicalGroup(1, [12], 7)
GridapGmsh.gmsh.model.setPhysicalName(1, 4, "bottom")
GridapGmsh.gmsh.model.setPhysicalName(1, 5, "right")
GridapGmsh.gmsh.model.setPhysicalName(1, 6, "top")
GridapGmsh.gmsh.model.setPhysicalName(1, 7, "left")



if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

gmsh.finalize()
