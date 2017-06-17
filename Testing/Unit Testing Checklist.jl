a1 = rand(4,6)
a2 = rand(5,6)
a3 = rand(3,6)
a4 = rand(5,6)
b1 = vcat(a1,a2,a3,a4)


function bar(x::Matrix, y...)
	vcat(x, y...)
end

c1 = Array(Matrix, 3)
c1[1] = a2
c1[2] = a3
c1[3] = a4
bar(a1, c1)

type MyType
	field1
end

mt = MyType(5)

## Unit Testing Checklist
#	- Population, Population - FitnessValues, and all other components
#		- concentrate on combine, duplicate and [], epecially multi-population combine
#	- EA_Parameters
#		- especially pμ_chr, pμ_xchr and pμ_nxchr
#	- EA_Counts
#		- especially select.reproduction and select.replication and the entire base/select fields
#		- updated partition counts (and slight changes to xover counts)
#