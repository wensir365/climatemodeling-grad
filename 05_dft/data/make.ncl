begin
	f 	= addfile("hgt.2017.nc","r")
	z	= f->hgt(332-1,{500},{55},:)
	write_matrix (onedtond(z,(/18,8/)), "8f10.2,", False)
end
