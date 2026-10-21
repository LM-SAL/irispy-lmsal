Removed ``irispy.utils.record_to_dict``; for a dictionary of a structured array's fields, use ``{name: arr[name] for name in arr.dtype.names}``.
