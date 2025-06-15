import numpy as np


class NeuralNetwork:
    def __init__(self, hidden_layers=(64, 64, 64), dim_input=1, learning_rate=0.001, epochs=5000):
        """        
        Parameters:
        - hidden_layers: tuple specifying number of neurons in each hidden layer
        - learning_rate: learning rate for gradient descent
        - epochs: number of training iterations
        """
        self.hidden_layers = hidden_layers
        self.learning_rate = learning_rate
        self.epochs = epochs
        self.weights = []
        self.biases = []
        self.dim = dim_input

    def initialise_parameters(self, input_dim):
        """Initialise weights and biases randomly."""
        layer_dims = [input_dim] + list(self.hidden_layers) + [self.dim]

        for i in range(len(layer_dims)-1):
            # Xavier/Glorot initialization
            bound = np.sqrt(6. / (layer_dims[i] + layer_dims[i+1]))
            self.weights.append(np.random.uniform(-bound, bound, (layer_dims[i], layer_dims[i+1])))
            self.biases.append(np.zeros((1, layer_dims[i+1])))

    def leaky_relu(self, x, alpha=0.01):
        """Leaky ReLU activation function to avoid dying neurons."""
        return np.where(x > 0, x, alpha * x)

    def leaky_relu_derivative(self, x, alpha=0.01):
        """Derivative of Leaky ReLU."""
        return np.where(x > 0, 1, alpha)

    def forward(self, X):
        """Forward propagation through the network."""
        activations = [X]
        pre_activations = []
    
        for i in range(len(self.weights)-1):
            z = np.dot(activations[-1], self.weights[i]) + self.biases[i]
            pre_activations.append(z)
            a = self.leaky_relu(z)
            activations.append(a)

        # Output layer (linear activation for regression)
        z = np.dot(activations[-1], self.weights[-1]) + self.biases[-1]
        pre_activations.append(z)
        activations.append(z)

        return activations, pre_activations

    def compute_loss(self, y_true, y_pred):
        """Mean squared error loss."""
        return np.mean((y_true - y_pred)**2)

    def backward(self, X, y, activations, pre_activations):
        """Backward propagation with gradient clipping."""
        m = X.shape[0]
        grads_w = [np.zeros_like(w) for w in self.weights]
        grads_b = [np.zeros_like(b) for b in self.biases]

        # Output layer gradient
        error = activations[-1] - y
        grads_w[-1] = np.dot(activations[-2].T, error) / m
        grads_b[-1] = np.sum(error, axis=0, keepdims=True) / m

        # Hidden layers gradients
        for l in range(len(self.weights)-2, -1, -1):
            error = np.dot(error, self.weights[l+1].T) * self.leaky_relu_derivative(pre_activations[l])
            # Gradient clipping to prevent exploding gradients
            error = np.clip(error, -1, 1)
            grads_w[l] = np.dot(activations[l].T, error) / m
            grads_b[l] = np.sum(error, axis=0, keepdims=True) / m

        return grads_w, grads_b

    def update_parameters(self, grads_w, grads_b):
        """Update weights and biases using gradient descent."""
        for i in range(len(self.weights)):
            self.weights[i] -= self.learning_rate * grads_w[i]
            self.biases[i] -= self.learning_rate * grads_b[i]

    def fit(self, X, y):
        """Train the neural network with learning rate scheduling."""
        n_samples, input_dim = X.shape
        self.initialise_parameters(input_dim)

        losses = []
        best_loss = float('inf')
        patience = 500
        patience_counter = 0

        for epoch in range(self.epochs):
            # Forward pass
            activations, pre_activations = self.forward(X)
            loss = self.compute_loss(y, activations[-1])
            losses.append(loss)

            # Early stopping check
            if loss < best_loss:
                best_loss = loss
                patience_counter = 0
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    print(f"Early stopping at epoch {epoch}")
                    break

            # Backward pass
            grads_w, grads_b = self.backward(X, y, activations, pre_activations)

            # Update parameters
            self.update_parameters(grads_w, grads_b)

            if epoch % 500 == 0:
                print(f"Epoch {epoch}, Loss: {loss:.6f}")

        return losses

    def predict(self, X):
        """Make predictions using the trained network."""
        activations, _ = self.forward(X)
        return activations[-1]
