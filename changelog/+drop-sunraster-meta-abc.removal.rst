``irispy.meta.SJIMeta`` and ``irispy.meta.SGMeta`` no longer have ``observer_location``, ``rsun_angular``, ``rsun_meters`` or ``version``, which were always ``None``.
For the Sun as seen from IRIS, use ``sun_angular_radius`` and ``distance_to_sun``; for the observer of a raster, ``SGMeta.observer``; for the data version, ``reformat_version`` or ``build_version``.
