package org.bigbio.pgatk.pepgenome.common;

public class SparkConfig {

    private static SparkConfig singleInstance = null;
    private String master = null;

    private SparkConfig() {
    }

    public static SparkConfig getInstance() {
        if (singleInstance == null) {
            singleInstance = new SparkConfig();
        }
        return singleInstance;
    }

    public String getMaster() {
        return master;
    }

    public void setMaster(String master) {
        this.master = master;
    }
}
