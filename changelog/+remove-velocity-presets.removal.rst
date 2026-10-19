Removed ``irispy.utils.rgb.linear_velocity`` and ``irispy.utils.rgb.asinh_velocity``; ``velocity_norm`` of `irispy.utils.rgb.calculate_rgb` and `irispy.utils.rgb.plot_rgb` now defaults to `None`, a linear mapping of the velocity that renders the same image.
For the arcsinh scale, pass ``velocity_norm=lambda v: np.arcsinh(v.to_value(u.km / u.s) / 25)``.
