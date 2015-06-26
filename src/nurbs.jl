using Gadfly
using Color
include("bspline.jl")
using bspline
x = linspace(0,1,1000)
Ξ = Float64[0.,0.,0.,0.5,1.,1.,1.]
H = Float64[0.,0.,0.,1.,1.,1.]
n1,p1 = bspline.knot_parameters(Ξ)
n2,p2 = bspline.knot_parameters(H)
N1(i,p) = Float64[ N(i,p,Ξ,x1)  for x1 in x ]
N2(i,p) = Float64[ N(i,p,H,x1)  for x1 in x ]


colors = distinguishable_colors(30)

plot(#[layer(x = x, y = N(i,1), Geom.line, Theme(default_color = colors[i])) for i in 1:8]...,
     [layer(x = x, y = N1(j,order), Geom.line, Theme(default_color = colors[j+2])) for j in 1:n1-p1+1]...,
     #[layer(x = x, y = N(k,3), Geom.line, Theme(default_color = colors[k+13])) for k in 1:5]...
     )

plot([layer(x = x, y = N2(j,order), Geom.line, Theme(default_color = colors[j+5])) for j in 1:n2-p2+1]...)

