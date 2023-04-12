def convert_units(value, units, desired_units=0):
    if units == "cb/a":
        value_in_desired_units = value * 0.148185
    elif units == "rb/a":
        value_in_desired_units = value * 13.6057
    elif units == "rb/cb":
        value_in_desired_units = value / 29421.02648438959
    else:
        value_in_desired_units = value
    return value_in_desired_units
