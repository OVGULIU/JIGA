using Gadfly
using Color
include("bspline.jl")
using bspline
x = linspace(0,1,1000)
y = linspace(0,1,1000)

Ξ = Float64[0.,0.,0.,0.5,1.,1.,1.]
H = Float64[0.,0.,0.,1.,1.,1.]

order = 2
n1 = bspline.knot_parameters(Ξ,order)
n2 = bspline.knot_parameters(H,order)

N1(i,p) = Float64[ bspline.N(i,p,Ξ,x1)  for x1 in x ]
M1(j,p) = Float64[ bspline.N(j,p,H,x1)  for x1 in x ]

bspline.global_basis_number(Ξ,H,1,3,order)
bspline.element_number(Ξ,H,n1,n2,2)

INC_array = bspline.INC(Ξ, H, order)
IEN_array = bspline.IEN(Ξ, H, order)

colors = distinguishable_colors(30)

plot(#[layer(x = x, y = N(i,1), Geom.line, Theme(default_color = colors[i])) for i in 1:8]...,
     [layer(x = x, y = N1(j,2), Geom.line, Theme(default_color = colors[j+2])) for j in 1:n1]...,
     #[layer(x = x, y = N(k,3), Geom.line, Theme(default_color = colors[k+13])) for k in 1:5]...
     )

plot(layer(x = x, y = N1(3,1), Geom.line, Theme(default_color = colors[5])))

