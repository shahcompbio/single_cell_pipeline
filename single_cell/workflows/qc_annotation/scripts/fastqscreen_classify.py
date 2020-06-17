import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import *


def train(training_data_path):
    '''
    Train the model using the provided training data.
    Return a feature scaler and a classifier.
    '''
    data = pd.read_csv(training_data_path)
    species_list = ["salmon", "grch37", "mm10"]
    labels = data["species"]
    features = data.drop('species', axis=1)

    le = LabelEncoder()
    le.fit(species_list)
    # convert the labels
    labels = le.transform(labels)
    # train a feature scaler
    transformer = RobustScaler().fit(features)
    features = transformer.transform(features)
    # train the random forest model
    rf = RandomForestClassifier(n_estimators=10, random_state=42)
    rf.fit(features, labels)

    return transformer, rf


def classify_fastqscreen(training_data_path, metrics_path, metrics_output):
    df = pd.read_csv(metrics_path)
    feature_transformer, model = train(training_data_path)

    features = ["fastqscreen_nohit_ratio", "fastqscreen_grch37_ratio", "fastqscreen_mm10_ratio",
                "fastqscreen_salmon_ratio"]
    label_to_species = {0: "grch37", 1: "mm10", 2: "salmon"}
    # check if all the features exists, if yes, make predictions, else create an empty species column.
    exist = all([feature[:-6] in df for feature in features])
    if exist:
        # make the feature columns
        for feature in features:
            df[feature] = df[feature[:-6]].divide(df["total_reads"])
        # scale the features
        scaled_features = feature_transformer.transform(df[features])
        df["species"] = model.predict(scaled_features)
        df["species"].replace(label_to_species, inplace=True)

    df.to_csv(metrics_output, index=False)
