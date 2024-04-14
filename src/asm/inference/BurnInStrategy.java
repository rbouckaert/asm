package asm.inference;

public enum BurnInStrategy {
    Running10("Running10"),
    Automatic("Automatic");

    private final String description;

    BurnInStrategy(String description) {
        this.description = description;
    }
//    private Integer burnInValue;

//    BurnInStrategy(String description, Integer value) {
//        this.description = description;
//        this.burnInValue = value;
//    }

//    public BurnInStrategy setBurnInValue(int newValue) {
//        this.burnInValue = newValue;
//        return this;
//    }

    public String getDescription() {
        return description;
    }

//    public Integer getBurnIn() {
//        return burnInValue;
//    }

    public static BurnInStrategy fromString(String input) {
//        String[] parts = input.split("-");
//        if (parts.length == 2) {
//            String name = parts[0];
//            int value = Integer.parseInt(parts[1]);
//            if (value >= 0 && value <= 99) {
//                for (BurnInStrategy enumConstant : BurnInStrategy.values()) {
//                    if (enumConstant.name().equalsIgnoreCase(name)) {
//                        // return new BurnInStrategy(name, value);
//                        BurnInStrategy newenum = enumConstant.valueOf(enumConstant.name());
//                        newenum.setBurnInValue(value);
//                        return newenum;
//                    }
//                }
////                return new ExampleEnum(name, value); // Create new enum with name and value
//            } else {
//                throw new IllegalArgumentException("Value must be between 0 and 99");
//            }
//        } else {
//            for (BurnInStrategy enumConstant : BurnInStrategy.values()) {
//                if (enumConstant.name().equalsIgnoreCase(input)) {
//                    return enumConstant;
//                }
//            }
//        }
//        throw new IllegalArgumentException("Invalid input: " + input);
        for (BurnInStrategy enumConstant : BurnInStrategy.values()) {
            if (enumConstant.name().equalsIgnoreCase(input)) {
                return enumConstant;
            }
        }
        throw new IllegalArgumentException("Invalid input: " + input);
    }




}
