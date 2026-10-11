Attribute VB_Name = "carvalho2025"
Option Explicit

Public Function LWAVE(ByVal PERIOD As Double, _
                      ByVal DEPTH As Double) As Variant
    Const g As Double = 9.80665
    Const PI As Double = 3.14159265358979
    Dim L0 As Double
    Dim alpha As Double
    Dim x As Double

    On Error GoTo InvalidInput
    If PERIOD <= 0# Or DEPTH <= 0# Then GoTo InvalidInput

    L0 = g * PERIOD ^ 2 / (2# * PI)
    alpha = 4# * PI ^ 2 * DEPTH / (g * PERIOD ^ 2)

    If alpha >= 20# Then
        LWAVE = L0
    Else
        x = 1.199315 ^ (alpha ^ 1.047086) _
            * alpha ^ 0.499947
        LWAVE = L0 * Application.WorksheetFunction.Tanh(x)
    End If
    Exit Function

InvalidInput:
    LWAVE = CVErr(xlErrNum)
End Function

