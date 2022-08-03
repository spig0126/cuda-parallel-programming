Errors
-------------
- Comparison between two floating-point numbers with the relational operator (==) is not accurate. Due to the internal precision errors in rounding upt floating-point nubmers, even though they have the same value, the comparison will say the two are not the same.