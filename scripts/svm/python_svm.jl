include("data_def.jl")
using ScikitLearn
@sk_import svm: SVC
model = SVC(kernel="rbf")
fit!(model, X_train, y_train);
accuracy = score(model, X_test, y_test)
y_pred_sk = predict(model, X_test);

function bench_python_svm(X_train, y_train)
  model = SVC(kernel="rbf")
  fit!(model, X_train, y_train);
end
