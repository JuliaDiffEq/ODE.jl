###############################################################################
## Coefficient Tableaus
###############################################################################

# Butcher Tableaus, or more generally coefficient tables
# see Hairer & Wanner 1992, p. 134, 166

abstract Tableau{T<:Real}
# Name is the name of the tableau/method (a symbol)
# S is the number of stages (an int)
# T is the type of the coefficients
#
# For all types of tableaus it assumes fields:
# order::(Int...) # order of the method(s)
#
# For Runge-Kutta methods it assumes fields:
# a::Matrix{T}  # SxS matrix
# b::Matrix{T}  # 1 or 2 x S matrix (fixed step/ adaptive)
# c::Vector{T}  # S
#
# For a tableau:
#  c1  | a_11   ....   a_1s
#  .   | a_21 .          .
#  .   | a_31     .      .
#  .   | ....         .  .
#  c_s | a_s1  ....... a_ss
# -----+--------------------
#      | b_1     ...   b_s   this is the one used for stepping
#      | b'_1    ...   b'_s  this is the one used for error-checking

Base.eltype{T}(b::Tableau{T}) = T
order(b::Tableau) = b.order
# Subtypes need to define a convert method to convert to a different
# eltype with signature:
Base.convert{T<:Tableau}(::Type{T}, tab::Tableau) = error("Define convert method for concrete Tableau types")
