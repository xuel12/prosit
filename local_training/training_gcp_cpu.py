import io_local
import losses
import model as model_lib
import constants
import constants_gcp_cpu

import os
from keras.utils import plot_model
import matplotlib.pyplot as plt
import numpy as np
import keras
from contextlib import redirect_stdout
import time

class TrainingPlot(keras.callbacks.Callback):
    def __init__(self, result_dir):
        self.result_dir = result_dir
        
    # This function is called when the training begins
    def on_train_begin(self, logs={}):
        # Initialize the lists for holding the logs, losses and accuracies
        self.losses = []
        self.acc = []
        self.val_losses = []
        self.val_acc = []
        self.logs = []

    # This function is called at the end of each epoch
    def on_epoch_end(self, epoch, logs={}):

        # Append the logs, losses and accuracies to the lists
        self.logs.append(logs)
        self.losses.append(logs.get('loss'))
        self.acc.append(logs.get('acc'))
        self.val_losses.append(logs.get('val_loss'))
        self.val_acc.append(logs.get('val_acc'))

        # Before plotting ensure at least 2 epochs have passed
        if len(self.losses) > 1:
            N = np.arange(0, len(self.losses))
            # You can chose the style of your preference
            # print(plt.style.available) to see the available options
            #plt.style.use("seaborn")

            # Plot train loss, train acc, val loss and val acc against epochs passed
            plt.figure()
            plt.plot(N, self.losses, label = "train_loss")
            plt.plot(N, self.acc, label = "train_acc")
            plt.plot(N, self.val_losses, label = "val_loss")
            plt.plot(N, self.val_acc, label = "val_acc")
            plt.title("Training Loss and Accuracy [Epoch {}]".format(epoch))
            plt.xlabel("Epoch #")
            plt.ylabel("Loss/Accuracy")
            plt.legend()
            # Make sure there exists a folder called output in the current directory
            # or replace 'output' with whatever direcory you want to put in the plots
            plt.savefig('{}/Epoch-{}.png'.format(self.result_dir, epoch))
            plt.close()


def get_callbacks(model_dir_path, result_dir):
    import keras

    loss_format = "{val_loss:.5f}"
    epoch_format = "{epoch:02d}"
    weights_file = "{}/weight_{}_{}.hdf5".format(
        model_dir_path, epoch_format, loss_format
    )
    csvlog_file = "{}/training.log".format(result_dir)
#    tensorboard = keras.callbacks.TensorBoard(log_dir='{}/tensorboardlogs'.format(result_dir), histogram_freq=1)
    save = keras.callbacks.ModelCheckpoint(weights_file, save_best_only=True)
    stop = keras.callbacks.EarlyStopping(patience=10)
    decay = keras.callbacks.ReduceLROnPlateau(patience=2, factor=0.2)
    csv_logger = keras.callbacks.CSVLogger(csvlog_file, append=False)
    plot_losses = TrainingPlot(result_dir)

#    return [save, stop, decay, csv_logger, plot_losses, tensorboard]
    return [save, stop, decay, csv_logger, plot_losses]


def train(tensor, model, model_config, callbacks):
#    import keras

    if isinstance(model_config["loss"], list):
        loss = [losses.get(l) for l in model_config["loss"]]
    else:
        loss = losses.get(model_config["loss"])
    optimizer = model_config["optimizer"]
    x = io_local.get_array(tensor, model_config["x"])
    y = io_local.get_array(tensor, model_config["y"])
    model.compile(optimizer=optimizer, loss=loss, metrics=['acc'])
    history = model.fit(
        x=x,
        y=y,
        epochs=constants.TRAIN_EPOCHS,
        batch_size=constants.TRAIN_BATCH_SIZE,
        validation_split=1 - constants.VAL_SPLIT,
        callbacks=callbacks,
    )
    keras.backend.get_session().close()
    return (history)
    

        
if __name__ == "__main__":
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"  # turn off tf logging
    os.chdir(constants_gcp_cpu.BASE_PATH + 'project/prosit/local_training')
    data_path = constants_gcp_cpu.DATA_PATH
    model_dir = constants_gcp_cpu.MODEL_DIR
    model, model_config = model_lib.load(model_dir, trained=False)
   
    # create log folder
    currenttime = int(time.time())
    result_dir = model_dir+'log_{}'.format(currenttime)   
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    # visualize model architecture
    plot_model(model, to_file='{}/model.png'.format(result_dir), show_shapes=True)
    with open('{}/modelsummary.txt'.format(result_dir), 'w') as f:
        with redirect_stdout(f):
            model.summary()
        
    tensor = io_local.from_hdf5(data_path)
    callbacks = get_callbacks(model_dir, result_dir)
    history = train(tensor, model, model_config, callbacks)
    history.history

    fig = plt.figure()  
    plt.plot(history.history['acc'])
    plt.plot(history.history['val_acc'])
    plt.title('Model accuracy')
    plt.ylabel('Accuracy')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Test'], loc='upper left')
#    plt.show()
    fig.savefig('{}/accuracy.png'.format(result_dir))
    
    # Plot training & validation loss values
    fig = plt.figure()  
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('Model loss')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Test'], loc='upper left')
#    plt.show()
    fig.savefig('{}/loss.png'.format(result_dir))

