subroutine help
print*,'This program is used to sub-box form icosahedral raw particle image'
print*,'icos_loc_box ...'
print*,'[icos_ort=ort0.dat]: icosahedral orientation and center file'
print*,'[local_ort=newort.dat]: output new orientation and center file'
print*,'<local_imgstck=substck.mrc>:output the stck file'
print*,'      if this parameter is not be defined, just output newort file'
print*,'<sub_FFTsize=256>'
print*,'<apix=1.0>'
print*,'<local_x=100.0> <local_y=100.0> <local_z=100.0>: (pixel)'
print*,'<PR_threshold=90>'
print*,'<boundX=100>'
print*,'<icos_imgstck=proj.stck>'
print*,'<sym_matrix=sym_maxtrix.txt>:symmetry operation matrix'
print*,'<first=1> <last=100>'
print*,'<sym=i2>: i2,i5,i3,c1,c2,...,c15'

print*,'OK, 2019-05-07. Hongrong Liu, HUNNU'

end
