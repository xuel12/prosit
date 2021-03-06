import numpy
import yaml
import sys
import os
try: 
    os.chdir('/Users/xuel12/Documents/MSdatascience/CS7180AI/project/prosit/local_training')
    print("Current directory is {}".format(os.getcwd()))
except: 
    print("Something wrong with specified directory. Exception- ", sys.exc_info())

import constants
import layers
import utils
from attention import Attention


def is_weight_name(w):
    return w.startswith("weight_") and w.endswith(".hdf5")


def get_loss(x):
    return float(x.split("_")[-1][:-5])


def get_best_weights_path(model_dir):
    weights = list(filter(is_weight_name, os.listdir(model_dir)))
    if len(weights) == 0:
        return None
    else:
        d = {get_loss(w): w for w in weights}
        weights_path = "{}/{}".format(model_dir, d[min(d)])
        return weights_path


def load(model_dir, trained=False):
    import keras

    model_path = os.path.join(model_dir, MODEL_NAME)
    config_path = os.path.join(model_dir, CONFIG_NAME)
    weights_path = get_best_weights_path(model_dir)
    with open(config_path, "r") as f:
        config = yaml.load(f.read(), Loader=yaml.SafeLoader)
    with open(model_path, "r") as f:
        model = keras.models.model_from_yaml(
            f.read(), custom_objects={"Attention": layers.Attention}
        )
    if trained and weights_path is not None:
        model.load_weights(weights_path)
    return model, config


def save(model, config, model_dir):
    model_path = os.path.join(model_dir, MODEL_NAME)
    config_path = os.path.join(model_dir, CONFIG_NAME)
    utils.check_mandatory_keys(config, ["name", "optimizer", "loss", "x", "y"])
    with open(config_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False)
    with open(model_path, "w") as f:
        f.write(model.to_yaml())
        

def model_build_mlp(modelfile, weightfile):
    from keras.models import Sequential 
    from keras.layers import Dense

    # fix random seed for reproducibility
    seed = 100
    numpy.random.seed(seed)
    
    # create model
    model = Sequential()
    model.add(Dense(12, input_dim=8, activation='relu'))
    model.add(Dense(8, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))
    
    # Compile model
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
     
    return model



def model_build_biGRU(modelfile, weightfile):
    from keras.models import Model
    from keras.layers import Input, LeakyReLU, Flatten, Dense, Dropout
    from keras.layers import Concatenate, Embedding, GRU, Bidirectional
    from keras.layers import RepeatVector, TimeDistributed, Multiply, Permute

    # fix random seed for reproducibility
    seed = 100
    numpy.random.seed(seed)
    
    peplen = 30
    max_features = 22
        
    # this embedding layer will encode the input sequence into a sequence of dense 32-dimensional vectors.
    peptides_in = Input(shape=(peplen,), dtype='int32', name='peptides_in')
    embedding = Embedding(max_features, 32, name='embedding')(peptides_in)
    encoder1 = Bidirectional(GRU(256, return_sequences=True, name = 'encoder1_gru'), name='encoder1')(embedding)
    dropout_1 = Dropout(0.3, name = 'dropout_1')(encoder1)
    encoder2 = GRU(512, return_sequences=True, name = 'encoder2')(dropout_1)
    dropout_2 = Dropout(0.3, name = 'dropout_2')(encoder2)
    encoder_att = Attention(name='encoder_att')(dropout_2)

    collision_energy_in = Input(shape=(1,), dtype='float32', name='collision_energy_in')
    precursor_charge_in = Input(shape=(6,), dtype='float32', name='precursor_charge_in')
    meta_in = Concatenate(axis=-1, name='meta_in')([collision_energy_in, precursor_charge_in])
    meta_dense = Dense(512, name='meta_dense')(meta_in)
    meta_dense_do = Dropout(0.3, name = 'meta_dense_do')(meta_dense)

    # combine seq, charge, ce embedding
    add_meta = Multiply(name='add_meta')([encoder_att, meta_dense_do])
    repeat = RepeatVector(29, name='repeat')(add_meta)
    decoder = GRU(512, return_sequences=True, name = 'decoder')(repeat)
    dropout_3 = Dropout(0.3, name = 'dropout_3')(decoder)
    
    permute_1 = Permute((2, 1), name = 'permute_1')(dropout_3)
    dense_1 = Dense(29, activation='softmax', name='dense_1')(permute_1)
    permute_2 = Permute((2, 1), name = 'permute_2')(dense_1)
    
    multiply_1 = Multiply(name='multiply_1')([dropout_3, permute_2])
    timedense = TimeDistributed(Dense(6, name='dense_2'), name='timedense')(multiply_1)
    activation = LeakyReLU(alpha=0.3, name = 'activation')(timedense) # names are added here
    out = Flatten(name='out')(activation)

    model = Model(input=[peptides_in, precursor_charge_in, collision_energy_in], output=[out], name='model_1')
    model.summary()
    
    return model

def model_build_biGRU_BN(modelfile, weightfile):
    from keras.models import Model
    from keras.layers import Input, LeakyReLU, Flatten, Dense, Dropout
    from keras.layers import Concatenate, Embedding, GRU, Bidirectional
    from keras.layers import RepeatVector, TimeDistributed, Multiply, Permute
    from keras.layers import BatchNormalization
    # fix random seed for reproducibility
    seed = 100
    numpy.random.seed(seed)
    
    peplen = 30
    max_features = 22
        
    # this embedding layer will encode the input sequence into a sequence of dense 32-dimensional vectors.
    peptides_in = Input(shape=(peplen,), dtype='int32', name='peptides_in')
    embedding = Embedding(max_features, 32, name='embedding')(peptides_in)
    encoder1 = Bidirectional(GRU(256, return_sequences=True, name = 'encoder1_gru'), name='encoder1')(embedding)
    encoder1 = BatchNormalization()(encoder1)
    dropout_1 = Dropout(0.3, name = 'dropout_1')(encoder1)
    encoder2 = GRU(512, return_sequences=True, name = 'encoder2')(dropout_1)
    dropout_2 = Dropout(0.3, name = 'dropout_2')(encoder2)
    encoder_att = Attention(name='encoder_att')(dropout_2)

    collision_energy_in = Input(shape=(1,), dtype='float32', name='collision_energy_in')
    precursor_charge_in = Input(shape=(6,), dtype='float32', name='precursor_charge_in')
    meta_in = Concatenate(axis=-1, name='meta_in')([collision_energy_in, precursor_charge_in])
    meta_dense = Dense(512, name='meta_dense')(meta_in)
    meta_dense_do = Dropout(0.3, name = 'meta_dense_do')(meta_dense)

    # combine seq, charge, ce embedding
    add_meta = Multiply(name='add_meta')([encoder_att, meta_dense_do])
    repeat = RepeatVector(29, name='repeat')(add_meta)
    decoder = GRU(512, return_sequences=True, name = 'decoder')(repeat)
    dropout_3 = Dropout(0.3, name = 'dropout_3')(decoder)
    
    permute_1 = Permute((2, 1), name = 'permute_1')(dropout_3)
    dense_1 = Dense(29, activation='softmax', name='dense_1')(permute_1)
    permute_2 = Permute((2, 1), name = 'permute_2')(dense_1)
    
    multiply_1 = Multiply(name='multiply_1')([dropout_3, permute_2])
    timedense = TimeDistributed(Dense(6, name='dense_2'), name='timedense')(multiply_1)
    activation = LeakyReLU(alpha=0.3, name = 'activation')(timedense) # names are added here
    out = Flatten(name='out')(activation)

    model = Model(input=[peptides_in, precursor_charge_in, collision_energy_in], output=[out], name='model_1')
    model.summary()
    
    return model

def model_build_biGRU_allBN(modelfile, weightfile):
    from keras.models import Model
    from keras.layers import Input, LeakyReLU, Flatten, Dense, Dropout
    from keras.layers import Concatenate, Embedding, GRU, Bidirectional
    from keras.layers import RepeatVector, TimeDistributed, Multiply, Permute
    from keras.layers import BatchNormalization
    # fix random seed for reproducibility
    seed = 100
    numpy.random.seed(seed)
    
    peplen = 30
    max_features = 22
        
    # this embedding layer will encode the input sequence into a sequence of dense 32-dimensional vectors.
    peptides_in = Input(shape=(peplen,), dtype='int32', name='peptides_in')
    embedding = Embedding(max_features, 32, name='embedding')(peptides_in)
    encoder1 = Bidirectional(GRU(256, return_sequences=True, name = 'encoder1_gru'), name='encoder1')(embedding)
    dropout_1 = BatchNormalization(name = 'batchnorm1')(encoder1)
#    dropout_1 = Dropout(0.3, name = 'dropout_1')(encoder1)
    encoder2 = GRU(512, return_sequences=True, name = 'encoder2')(dropout_1)
    dropout_2 = BatchNormalization(name = 'batchnorm2')(encoder2)
#    dropout_2 = Dropout(0.3, name = 'dropout_2')(encoder2)
    encoder_att = Attention(name='encoder_att')(dropout_2)

    collision_energy_in = Input(shape=(1,), dtype='float32', name='collision_energy_in')
    precursor_charge_in = Input(shape=(6,), dtype='float32', name='precursor_charge_in')
    meta_in = Concatenate(axis=-1, name='meta_in')([collision_energy_in, precursor_charge_in])
    meta_dense = Dense(512, name='meta_dense')(meta_in)
    meta_dense_do = BatchNormalization(name = 'meta_dense_do')(meta_dense)
#    meta_dense_do = Dropout(0.3, name = 'meta_dense_do')(meta_dense)

    # combine seq, charge, ce embedding
    add_meta = Multiply(name='add_meta')([encoder_att, meta_dense_do])
    repeat = RepeatVector(29, name='repeat')(add_meta)
    decoder = GRU(512, return_sequences=True, name = 'decoder')(repeat)
    dropout_3 = BatchNormalization(name = 'batchnorm3')(decoder)
#    dropout_3 = Dropout(0.3, name = 'dropout_3')(decoder)
    
    permute_1 = Permute((2, 1), name = 'permute_1')(dropout_3)
    dense_1 = Dense(29, activation='softmax', name='dense_1')(permute_1)
    permute_2 = Permute((2, 1), name = 'permute_2')(dense_1)
    
    multiply_1 = Multiply(name='multiply_1')([dropout_3, permute_2])
    timedense = TimeDistributed(Dense(6, name='dense_2'), name='timedense')(multiply_1)
    timedense = BatchNormalization(name = 'batchnorm4')(timedense)
    activation = LeakyReLU(alpha=0.3, name = 'activation')(timedense) # names are added here
    out = Flatten(name='out')(activation)

    model = Model(input=[peptides_in, precursor_charge_in, collision_energy_in], output=[out], name='model_1')
    model.summary()
    
    return model


def model_build_mlp(modelfile, weightfile):
    from keras.models import Model, Sequential
    from keras.layers import Input, LeakyReLU, Flatten, Dense, Dropout
    from keras.layers import Concatenate, Embedding, GRU, Bidirectional
    from keras.layers import RepeatVector, TimeDistributed, Multiply, Permute
    from keras.layers import BatchNormalization
    # fix random seed for reproducibility
    seed = 100
    numpy.random.seed(seed)
    
    peplen = 30
#    max_features = 22
    
    # this embedding layer will encode the input sequence into a sequence of dense 32-dimensional vectors.
    peptides_in = Input(shape=(peplen,), dtype='float32', name='peptides_in')
    collision_energy_in = Input(shape=(1,), dtype='float32', name='collision_energy_in')
    precursor_charge_in = Input(shape=(6,), dtype='float32', name='precursor_charge_in')
    meta_in = Concatenate(axis=-1, name='meta_in')([peptides_in, collision_energy_in, precursor_charge_in])
    meta_dense = Dense(512, name='meta_dense')(meta_in)
    hidden1 = Dense(512,activation='relu', name = 'hidden1')(meta_dense)
    hidden2 = Dense(512,activation='relu', name = 'hidden2')(hidden1)
    dropout_2 = BatchNormalization(name = 'batchnorm1')(hidden2)
    hidden3 = Dense(512,activation='relu', name = 'hidden3')(dropout_2)
    hidden4 = Dense(512,activation='relu', name = 'hidden4')(hidden3)
    dropout_3 = Dropout(0.1, name = 'dropout_3')(hidden4)
    hidden5 = Dense(256,activation='relu',name = 'hidden5')(dropout_3)
    hidden6 = Dense(256,activation='relu',name = 'hidden6')(hidden5)
    hidden7 = Dense(256,activation='relu',name = 'hidden7')(hidden6)
    dropout_3 = BatchNormalization(name = 'batchnorm2')(hidden7)
    activation = LeakyReLU(alpha=0.3, name = 'activation')(dropout_3) # names are added here
    out = Dense(174, name='out')(activation)
    model = Model(input=[peptides_in, precursor_charge_in, collision_energy_in], output=[out], name='model_1')
    model.summary()
    
    return model


def model_build_seq2seq(modelfile, weightfile):
    from keras.models import Model
    from keras.layers import Input, LeakyReLU, Flatten, Dense, Dropout
    from keras.layers import Concatenate, Embedding, GRU, Bidirectional,SimpleRNN,LSTM
    from keras.layers import RepeatVector, TimeDistributed, Multiply, Permute

    # fix random seed for reproducibility
    seed = 100
    numpy.random.seed(seed)
    
    peplen = 30
    max_features = 22
        
    # this embedding layer will encode the input sequence into a sequence of dense 32-dimensional vectors.
    peptides_in = Input(shape=(peplen,), dtype='int32', name='peptides_in')
    embedding = Embedding(max_features, 32, name='embedding')(peptides_in)

    collision_energy_in = Input(shape=(1,), dtype='float32', name='collision_energy_in')
    precursor_charge_in = Input(shape=(6,), dtype='float32', name='precursor_charge_in')
    meta_in = Concatenate(axis=-1, name='meta_in')([collision_energy_in, precursor_charge_in])
    meta_dense = RepeatVector(30, name='meta_dense')(meta_in)
    meta_dense_do = Dense(32, name='meta_dense_do')(meta_dense)
    add_meta = Multiply(name='add_meta')([embedding, meta_dense_do])

#    encoder1 = Bidirectional(GRU(256, return_sequences=True, name = 'encoder1_gru'), name='encoder1')(add_meta)
    encoder1 = Bidirectional(LSTM(256, return_sequences=True, name = 'encoder1_lstm'), name='encoder1')(add_meta)
#    encoder1 = GRU(512, return_sequences=True, name = 'encoder1_gru')(add_meta)
#    encoder1 = LSTM(512, return_sequences=True, name = 'encoder1_lstm')(add_meta)
#    encoder1 = SimpleRNN(512, return_sequences=True, name = 'encoder1_rnn')(add_meta)

    dropout_1 = Dropout(0.3, name = 'dropout_1')(encoder1)
    encoder_att = Attention(name='encoder_att')(dropout_1)

    repeat = RepeatVector(29, name='repeat')(encoder_att)
#    decoder = Bidirectional(GRU(256, return_sequences=True, name = 'decoder_gru'), name='decoder')(repeat)
    decoder = Bidirectional(LSTM(256, return_sequences=True, name = 'decoder_lstm'), name='decoder')(repeat)
#    decoder = GRU(512, return_sequences=True, name = 'decoder')(repeat)
#    decoder = LSTM(512, return_sequences=True, name = 'decoder')(repeat)
#    decoder = SimpleRNN(512, return_sequences=True, name = 'decoder')(repeat)

#    
#    permute_1 = Permute((2, 1), name = 'permute_1')(dropout_3)
#    dense_1 = Dense(29, activation='softmax', name='dense_1')(permute_1)
#    permute_2 = Permute((2, 1), name = 'permute_2')(dense_1)
#    
#    multiply_1 = Multiply(name='multiply_1')([dropout_3, permute_2])
    timedense = Dense(6,activation='relu',name = 'timedense')(decoder)

#    timedense = TimeDistributed(Dense(6, name='dense_2'), name='timedense')(decoder)
    activation = LeakyReLU(alpha=0.3, name = 'activation')(timedense) # names are added here
    out = Flatten(name='out')(activation)

    model = Model(input=[peptides_in, precursor_charge_in, collision_energy_in], output=[out], name='model_1')
    model.summary()
    
    return model

def model_build_LSTM(modelfile, weightfile):
    from keras.models import Model
    from keras.layers import Input, LeakyReLU, Flatten, Dense, Dropout
    from keras.layers import Concatenate, Embedding, GRU, Bidirectional,SimpleRNN,LSTM,CuDNNLSTM
    from keras.layers import RepeatVector, TimeDistributed, Multiply, Permute

    # fix random seed for reproducibility
    seed = 100
    numpy.random.seed(seed)
    
    peplen = 30
    max_features = 22
        
    # this embedding layer will encode the input sequence into a sequence of dense 32-dimensional vectors.
    peptides_in = Input(shape=(peplen,), dtype='int32', name='peptides_in')
    embedding = Embedding(max_features, 32, name='embedding')(peptides_in)
    encoder1 = CuDNNLSTM(512, return_sequences=True, name = 'encoder1')(embedding)
#    dropout_1 = Dropout(0.3, name = 'dropout_1')(encoder1)
#    encoder2 = LSTM(512, return_sequences=True, name = 'encoder2')(dropout_1)
    dropout_2 = Dropout(0.3, name = 'dropout_2')(encoder1)
    encoder_att = Attention(name='encoder_att')(dropout_2)

    collision_energy_in = Input(shape=(1,), dtype='float32', name='collision_energy_in')
    precursor_charge_in = Input(shape=(6,), dtype='float32', name='precursor_charge_in')
    meta_in = Concatenate(axis=-1, name='meta_in')([collision_energy_in, precursor_charge_in])
    meta_dense = Dense(512, name='meta_dense')(meta_in)
    meta_dense_do = Dropout(0.3, name = 'meta_dense_do')(meta_dense)

    # combine seq, charge, ce embedding
    add_meta = Multiply(name='add_meta')([encoder_att, meta_dense_do])
    repeat = RepeatVector(29, name='repeat')(add_meta)
    decoder = CuDNNLSTM(512, return_sequences=True, name = 'decoder')(repeat)
    dropout_3 = Dropout(0.3, name = 'dropout_3')(decoder)
    
    permute_1 = Permute((2, 1), name = 'permute_1')(dropout_3)
    dense_1 = Dense(29, activation='softmax', name='dense_1')(permute_1)
    permute_2 = Permute((2, 1), name = 'permute_2')(dense_1)
    
    multiply_1 = Multiply(name='multiply_1')([dropout_3, permute_2])
    timedense = TimeDistributed(Dense(6, name='dense_2'), name='timedense')(multiply_1)
    activation = LeakyReLU(alpha=0.3, name = 'activation')(timedense) # names are added here
    out = Flatten(name='out')(activation)

    model = Model(input=[peptides_in, precursor_charge_in, collision_energy_in], output=[out], name='model_1')
    model.summary()
    
    return model


if __name__ == "__main__":
    
    os.chdir(constants.BASE_PATH + 'project/prosit/local_training')
    model_dir = constants.MODEL_DIR
    MODEL_NAME = "model.yml"
    CONFIG_NAME = "config.yml"
    WEIGHT_NAME = "weight.hdf5"
    
    model, config = load(model_dir, trained=False)
#    model = model_build_biGRU(model_dir + MODEL_NAME, model_dir + WEIGHT_NAME)
    model = model_build_LSTM(model_dir + MODEL_NAME, model_dir + WEIGHT_NAME)

    save(model, config, model_dir+'/tmpmodel')
    
