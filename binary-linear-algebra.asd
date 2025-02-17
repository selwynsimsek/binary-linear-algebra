(defsystem "binary-linear-algebra"
  :version "0.0.1"
  :author "Selwyn Simsek"
  :license "BSD 2"
  :depends-on ("metabang-bind"
               "zr-utils"
               "iterate")
  :components ((:module "src"
                :components
                ((:file "package")
                 (:file "routines")
                 (:file "aliases"))))
  :description "Linear algebra over 𝔽₂ in Common Lisp"
  :in-order-to ((test-op (test-op "binary-linear-algebra/tests"))))

(defsystem "binary-linear-algebra/tests"
  :author "Selwyn Simsek"
  :license "BSD 2"
  :depends-on ("binary-linear-algebra"
               "rove")
  :components ((:module "tests"
                :components
                ((:file "main"))))
  :description "Test system for binary-linear-algebra"
  :perform (test-op (op c) (symbol-call :rove :run c)))
