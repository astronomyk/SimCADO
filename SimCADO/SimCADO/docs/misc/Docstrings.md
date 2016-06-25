# Quick and dirty guide to Docstrings

Taken from the Astropy Documentation guide

## Functions
1. One-line summary
1. Long summary - clarify functionality
1. Parameters
```
x  :  type [, optional [, {set values} ]]
	Description  of  parameter  `x`. [(Default value)]
```

4. Returns
1. Keyword Arguments (**kwargs)
1. Raises
```
InvalidWCSException
	If  the  WCS  information  is  invalid.
```

7. See Also - point to functions of classes that are connect
1. Notes
1. References
1. Examples


```
"""
<One-line summary goes here>


Summary
-------


Parameters
----------
x  :  type [, optional [, {set values} ]]
	Description of `x`. [(Default value)]

Returns
-------


Keyword Arguments (**kwargs)
----------------------------


Raises
------


See Also
--------


Notes
-----


References
----------


Examples
--------
"""
```

---


## Classes
1. One-line summary
1. Long summary - clarify functionality
1. Parameters
1. Attributes
1. Methods - Private methods go in Notes
1. Keyword Arguments (**kwargs)
1. Raises
1. See Also - point to functions of classes that are connect
1. Notes
1. References
1. Examples

```
"""
<One-line summary goes here>


Summary
-------


Parameters
----------


Attributes
----------


Methods
-------


Keyword Arguments (**kwargs)
----------------------------


Raises
------


See Also
--------


Notes
-----


References
----------


Examples
--------
"""
```


## Modules
1. One-line summary
1. Long summary
1. Routine  listings
1. See  Also
1. Notes
1. References
1. Examples


```
"""
<One-line summary goes here>


Summary
-------


Routines
--------


See Also
--------


Notes
-----


References
----------


Examples
--------
"""
```

