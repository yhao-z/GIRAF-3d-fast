# GIRAF-3d-fast
GIRAF algorithm for MRI data in 3 dimension without the CG iterative algorithm, using two steps of variable splitting instead of only one step in the [the original 2d GIRAF algorithm](https://github.com/cbig-iowa/giraf), that means the algorithm is fast.

the codes of this method is just a 3d extension of [the original 2d GIRAF algorithm](https://github.com/cbig-iowa/giraf) from [Greg Ongie](https://github.com/gregongie)

because of the 3d properties, it is difficult to using the same solution of the original optimization problem. thus, I use a two steps of variable splitting instead, resulting in an analytical form for all the iterative steps instead of importing CG methods.

the derivation of formulas will be upload as soon as possible.

for now, it is an early version, because the annotation of the codes is not so well and is mixed with English and Chinese. it will be correct in a short time.

although there are few deficiencies, it is easily to run the codes correctly. try it and test it!

perfect it together!

# NEXT IS THE derivation
download and look the pdf file [Derivation of GIRAF-3D-FAST.pdf](https://github.com/yhao-z/GIRAF-3d-fast/blob/main/Derivation%20of%20GIRAF-3D-FAST.pdf) for detail.
