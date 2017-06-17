import Base: <=, <, >=, >, ==

####################################################################################################
# <= comparisons
####################################################################################################

function <=(x::Vector, y::Vector)
	lngth = length(x)
	if lngth != length(y) error("Function <=(::Vector, ::Vector), vectors must be same length and are not") end
	ans = Array(Bool, lngth)
	for i = 1:lngth
		ans[i] = x[i] <= y[i]
	end
	ans
end

function <=(xv::Vector, y::Real)
	lngth = length(xv)
	ans = Array(Bool, lngth)
	yv = fill(y, lngth)
	xv <= yv
end

function <=(x::Real, yv::Vector)
	lngth = length(yv)
	ans = Array(Bool, lngth)
	xv = fill(x, lngth)
	xv <= yv
end

####################################################################################################
# < comparisons
####################################################################################################

function <(x::Vector, y::Vector)
	lngth = length(x)
	if lngth != length(y) error("Function <=(::Vector, ::Vector), vectors must be same length and are not") end
	ans = Array(Bool, lngth)
	for i = 1:lngth
		ans[i] = x[i] < y[i]
	end
	ans
end

function <(xv::Vector, y::Real)
	lngth = length(xv)
	ans = Array(Bool, lngth)
	yv = fill(y, lngth)
	xv < yv
end

function <(x::Real, yv::Vector)
	lngth = length(yv)
	ans = Array(Bool, lngth)
	xv = fill(x, lngth)
	xv < yv
end

####################################################################################################
# >= comparisons
####################################################################################################

function >=(x::Vector, y::Vector)
	lngth = length(x)
	if lngth != length(y) error("Function <=(::Vector, ::Vector), vectors must be same length and are not") end
	ans = Array(Bool, lngth)
	for i = 1:lngth
		ans[i] = x[i] >= y[i]
	end
	ans
end

function >=(xv::Vector, y::Real)
	lngth = length(xv)
	ans = Array(Bool, lngth)
	yv = fill(y, lngth)
	xv >= yv
end

function >=(x::Real, yv::Vector)
	lngth = length(yv)
	ans = Array(Bool, lngth)
	xv = fill(x, lngth)
	xv >= yv
end

####################################################################################################
# > comparisons
####################################################################################################

function >(x::Vector, y::Vector)
	lngth = length(x)
	if lngth != length(y) error("Function <=(::Vector, ::Vector), vectors must be same length and are not") end
	ans = Array(Bool, lngth)
	for i = 1:lngth
		ans[i] = x[i] > y[i]
	end
	ans
end

function >(xv::Vector, y::Real)
	lngth = length(xv)
	ans = Array(Bool, lngth)
	yv = fill(y, lngth)
	xv > yv
end

function >(x::Real, y::Vector)
	lngth = length(y)
	ans = Array(Bool, lngth)
	xv = fill(x, lngth)
	xv > y
end

####################################################################################################
# == comparisons
####################################################################################################

function ==(x::Vector, y::Vector)
	lngth = length(x)
	if lngth != length(y) error("Function <=(::Vector, ::Vector), vectors must be same length and are not") end
	ans = Array(Bool, lngth)
	for i = 1:lngth
		ans[i] = x[i] == y[i]
	end
	ans
end

function ==(xv::Vector, y::Number)
	lngth = length(xv)
	ans = Array(Bool, lngth)
	yv = fill(y, lngth)
	xv == yv
end

function ==(x::Number, yv::Vector)
	lngth = length(yv)
	ans = Array(Bool, lngth)
	xv = fill(x, lngth)
	xv == yv
end
