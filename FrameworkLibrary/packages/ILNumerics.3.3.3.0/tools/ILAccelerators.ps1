
Import-Module TypeAccelerator

Add-TypeAccelerator ilf32 ILNumerics.ILArray[single]
Add-TypeAccelerator ilif32 ILNumerics.ILInArray[single]
Add-TypeAccelerator ilof32 ILNumerics.ILOutArray[single]
Add-TypeAccelerator ilrf32 ILNumerics.ILRetArray[single]

Add-TypeAccelerator ilf64 ILNumerics.ILArray[double]
Add-TypeAccelerator ilif64 ILNumerics.ILInArray[double]
Add-TypeAccelerator ilof64 ILNumerics.ILOutArray[double]
Add-TypeAccelerator ilrf64 ILNumerics.ILRetArray[double]

Add-TypeAccelerator ili8 ILNumerics.ILArray[SByte]
Add-TypeAccelerator ilii8 ILNumerics.ILInArray[SByte]
Add-TypeAccelerator iloi8 ILNumerics.ILOutArray[SByte]
Add-TypeAccelerator ilri8 ILNumerics.ILRetArray[SByte]

Add-TypeAccelerator ili16 ILNumerics.ILArray[Int16]
Add-TypeAccelerator ilii16 ILNumerics.ILInArray[Int16]
Add-TypeAccelerator iloi16 ILNumerics.ILOutArray[Int16]
Add-TypeAccelerator ilri16 ILNumerics.ILRetArray[Int16]

Add-TypeAccelerator ili32 ILNumerics.ILArray[Int32]
Add-TypeAccelerator ilii32 ILNumerics.ILInArray[Int32]
Add-TypeAccelerator iloi32 ILNumerics.ILOutArray[Int32]
Add-TypeAccelerator ilri32 ILNumerics.ILRetArray[Int32]

Add-TypeAccelerator ilu8 ILNumerics.ILArray[Byte]
Add-TypeAccelerator iliu8 ILNumerics.ILInArray[Byte]
Add-TypeAccelerator ilou8 ILNumerics.ILOutArray[Byte]
Add-TypeAccelerator ilru8 ILNumerics.ILRetArray[Byte]

Add-TypeAccelerator ilstring ILNumerics.ILArray[String]
Add-TypeAccelerator ilistring ILNumerics.ILInArray[String]
Add-TypeAccelerator ilostring ILNumerics.ILOutArray[String]
Add-TypeAccelerator ilrstring ILNumerics.ILRetArray[String]

Add-TypeAccelerator ilcell ILNumerics.ILCell
Add-TypeAccelerator ilicell ILNumerics.ILInCell
Add-TypeAccelerator ilocell ILNumerics.ILOutCell
Add-TypeAccelerator ilrcell ILNumerics.ILRetCell

Add-TypeAccelerator ilbool ILNumerics.ILLogical
Add-TypeAccelerator ilibool ILNumerics.ILInLogical
Add-TypeAccelerator ilobool ILNumerics.ILOutLogical
Add-TypeAccelerator ilrbool ILNumerics.ILRetLogical

Add-TypeAccelerator ilmath ILNumerics.ILMath
Add-TypeAccelerator ilsettings ILNumerics.Settings
Add-TypeAccelerator ilscope ILNumerics.ILScope
