package asm.inference;

public enum BurnInStrategy {
    Running("Running"),
    Automatic("Automatic");

    private final String strategy;

    private static final int DEFAULT_BURN_IN_PERCENT = 10;

    public int getBurnInPercent() {
        return burnInPercent;
    }

    int burnInPercent;

    BurnInStrategy(String strategy) {
        this.strategy = strategy;
        if (strategy.equalsIgnoreCase("Automatic")) {
            this.burnInPercent = 0;
        } else if (strategy.equalsIgnoreCase("Running")) {
            this.burnInPercent = DEFAULT_BURN_IN_PERCENT;
        }
    }

    public String getStrategy() {
        return strategy;
    }

    void setBurnInPercent(int burnInPercent) {
        if (burnInPercent >= 0 && burnInPercent <= 99) {
            this.burnInPercent = burnInPercent;
        } else {
            throw new IllegalArgumentException("BurnInPercent must be between 0 and 99");
        }
    }
}
