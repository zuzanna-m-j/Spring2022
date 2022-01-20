/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 2.0.9
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package libgetar_wrap;

public final class Format {
  public final static Format Float32 = new Format("Float32");
  public final static Format Float64 = new Format("Float64");
  public final static Format Int32 = new Format("Int32");
  public final static Format Int64 = new Format("Int64");
  public final static Format UInt8 = new Format("UInt8");
  public final static Format UInt32 = new Format("UInt32");
  public final static Format UInt64 = new Format("UInt64");

  public final int swigValue() {
    return swigValue;
  }

  public String toString() {
    return swigName;
  }

  public static Format swigToEnum(int swigValue) {
    if (swigValue < swigValues.length && swigValue >= 0 && swigValues[swigValue].swigValue == swigValue)
      return swigValues[swigValue];
    for (int i = 0; i < swigValues.length; i++)
      if (swigValues[i].swigValue == swigValue)
        return swigValues[i];
    throw new IllegalArgumentException("No enum " + Format.class + " with value " + swigValue);
  }

  private Format(String swigName) {
    this.swigName = swigName;
    this.swigValue = swigNext++;
  }

  private Format(String swigName, int swigValue) {
    this.swigName = swigName;
    this.swigValue = swigValue;
    swigNext = swigValue+1;
  }

  private Format(String swigName, Format swigEnum) {
    this.swigName = swigName;
    this.swigValue = swigEnum.swigValue;
    swigNext = this.swigValue+1;
  }

  private static Format[] swigValues = { Float32, Float64, Int32, Int64, UInt8, UInt32, UInt64 };
  private static int swigNext = 0;
  private final int swigValue;
  private final String swigName;
}
