from proj_4 import Diffusion_1D, Diffusion_3D

bar = Diffusion_1D()
bar.set_boundary_conditions()
bar.simulate_diffusion()

sphere = Diffusion_3D(t_max = 180)
sphere.set_boundary_conditions()
sphere.simulate_diffusion()
sphere.convert_temp_to_radial()

bar.plot_animate("1D bar")
sphere.plot_animate("3D sphere")
