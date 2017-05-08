%mem=3gb
%chk=calculation_0002
%nprocshared=4
#p MP4 6-311+G(d) nosym field=read
scf=(Conver=11,NoVarAcc,MaxCyc=600,vshift=1000) IOP(9/6=600,9/9=11)

finite field calculation (-1, 0, 0)

-1 2
O	  0.00000000   0.00000000   0.00000000
H	  0.75698805   0.00000000   0.58632839
H	 -0.75698805   0.00000000   0.58632839

-0.0004000000	 0.0000000000	 0.0000000000