#Firzandi Ahsan DS-21537144016
function diff(x::Array,y::Array)
	m=length(x) 
	a=Array{Float64}(undef,m)
	for i in 1:m
		a[i]=y[i]
	end
	for j in 2:m
		for i in reverse(collect(j:m))
			a[i]=(a[i]-a[i-1])/(x[i]-x[i-(j-1)])
		end
	end
	return(a)
end

function newton(x::Array,y::Array,z)
	m=length(x) 
	a=diff(x,y)
	sum=a[1]
	pr=1.0
	for j in 1:(m-1)
		pr=pr*(z-x[j])
		sum=sum+a[j+1]*pr
	end
	return sum
end


#Firzandi Ahsan DS-21537144016
using PyPlot
using LaTeXStrings

f(x)=x*sin(x)
xi-collect(-5:10/3:5) 
yi=map(f,xi) 
xaxis=-5:1/100:5
runge=map(f, xaxis)
interp=map(z->newton (xi, yi,z), xaxis) 
plot (xaxis, interp, label="interpolating poly")
plot (xaxis, runge, label=L"f(x)=x*sin(x)")
scatter(xi, yi, label="data")
legend (loc="upper right");


#Firzandi Ahsan DS-21537144016
using Pyplot
using LaTeXStrings

f(x)=x*sin(x)
xi-collect(-pi: 10/3: pi)  
yi=map(f,xi) 
xaxis=-5:1/100:5 plot (xaxis, interp, label="interpolating poly")
runge=map(f, xaxis) 
interp=map(z->newton (xi, yi,z), xaxis) 
plot(xaxis, interp, label="interpolating poly")
plot(xaxis, runge, label=L"f(x)=x*sin(x)")
scatter(xi, yi, label="data")
legend (loc= "lower center");
savefig("coba1.png")


#Firzandi Ahsan DS-21537144016
xi-collect(-pi: 10/5:pi) # 6 equally spaced values from -5 to 5
yi=map(f,xi)
interp=map(z->newton (xi,yi,z), xaxis) 
plot (xaxis, interp, label="interpolating poly") 
plot (xaxis, runge,label=L"f(x)=1/(1+x^2)")
scatter (xi, yi, label="data") 
legend (loc= "lower center"); 
savefig("coba2.png")


#Firzandi Ahsan DS-21537144016
xi-collect(-pi: 10/10: pi) 
yi=map(f,xi) 
interp=map(z->newton (xi,yi,z), xaxis) 
plot(xaxis, interp, label="interpolating poly") 
plot(xaxis, runge, label=L"f(x)=x*sin(x)") 
scatter (xi, yi, label="data") 
legend (loc= "upper center");


#Firzandi Ahsan DS-21537144016
xi-collect(-pi: 10/20: pi) 
yi=map(f,xi)
interp=map (z->newton (xi, yi,z), xaxis) 
plot(xaxis, interp, label="interpolating poly") 
plot(xaxis, runge, label=L"f(x)=x*sin(x)") 
scatter(xi, yi, label="data") 
legend (loc= "lower center"); 
savefig("coba4.png")