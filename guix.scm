;; To use this file to build a version of moni using git HEAD:
;;
;;  rm -rf build/*
;;  guix build -f guix.scm
;;
;; To get a development container:
;;
;;   guix shell --share=/path/to/host/folder=/path/to/mount/inside/guix --network -C -D -f guix.scm -- bash --init-file <(echo "mkdir -p /usr/bin && ln -s \$GUIX_ENVIRONMENT/bin/env /usr/bin/env")
;;
;; and inside the container:
;;
;;   mkdir build
;;   cd build
;;   cmake ..
;;   make
;;
;; You can then run moni inside the container.
;;
;; Andrea Guarracino (c) 2023-2024

(use-modules
 ((guix licenses) #:prefix license:)
 (guix gexp)
 (guix packages)
 (guix git-download)
 (guix build-system cmake)
 (guix utils)
 (gnu packages base)
 (gnu packages compression)
 (gnu packages curl)
 (gnu packages gcc)
 (gnu packages python)
 (gnu packages version-control)
 (gnu packages tls)
 (gnu packages ncurses)
 (gnu packages autotools)
 (gnu packages cmake)
 (gnu packages pkg-config)
 (gnu packages boost)
 (gnu packages certs))

(define %source-dir (dirname (current-filename)))

(define-public moni-align-git
  (package
    (name "moni-align-git")
    (version "git")
    (source (origin
              (method git-fetch)
              (uri (git-reference
                    (url "https://github.com/maxrossi91/moni-align.git")
                    (commit "develop")))
              (file-name (git-file-name name version))
              (sha256
               (base32
                "0fx0ji7hk9hgskl3psm68f836vy7pxny37y2lq9a7pm456gz2wrd"))))
    (build-system cmake-build-system)
    (arguments
     `(#:tests? #f
       #:configure-flags '("-DCMAKE_BUILD_TYPE=Release"
                           "-DCMAKE_C_FLAGS=-fcommon"
                           "-DCMAKE_CXX_FLAGS=-fcommon")
       #:phases
       (modify-phases %standard-phases
         (add-after 'unpack 'fix-newscan-cpp
           (lambda _
             (let ((file "thirdparty/bigrepair/ctph/newscan.cpp"))
               (chmod file #o644)
               (substitute* file
                 (("stoi\\( sarg \\)")
                  "std::stoi(std::string(sarg))")))))
         (add-after 'unpack 'set-env
           (lambda* (#:key inputs outputs #:allow-other-keys)
             (when (getenv "LD_LIBRARY_PATH")
               (setenv "LD_LIBRARY_PATH"
                       (string-append (getenv "LD_LIBRARY_PATH")
                                      ":/build/moni-align/build/thirdparty/lib")))
             (setenv "CC" (which "gcc"))
             (setenv "CXX" (which "g++"))
             (setenv "SSL_CERT_FILE" (string-append (assoc-ref inputs "nss-certs")
                                                    "/etc/ssl/certs/ca-certificates.crt"))
             (setenv "PATH" (string-append (assoc-ref outputs "out") "/bin:"
                                           (getenv "PATH")))
             #t))
         (add-after 'unpack 'fix-cmake-modules-path
           (lambda* (#:key inputs #:allow-other-keys)
             (let ((cmake-root (assoc-ref inputs "cmake")))
               (when cmake-root
                 (setenv "CMAKE_PREFIX_PATH"
                         (string-append cmake-root "/share/cmake-"
                                        ,(package-version cmake) "/Modules"))))
             #t))
         (add-before 'configure 'fix-thirdparty-paths
           (lambda _
             (substitute* "CMakeLists.txt"
               (("\\$\\{CMAKE_SOURCE_DIR\\}/thirdparty")
                (string-append (getcwd) "/thirdparty")))
             #t))
         (add-after 'build 'build-thirdparty
           (lambda* (#:key outputs #:allow-other-keys)
             (let ((out (assoc-ref outputs "out")))
               (with-directory-excursion "thirdparty"
                 ;; Build SDSL
                 (invoke "cmake" "-DCMAKE_INSTALL_PREFIX=${out}" "sdsl-lite")
                 (invoke "make")
                 (invoke "make" "install")
                 
                 ;; Build htslib
                 (with-directory-excursion "htslib"
                   (invoke "autoreconf" "-i")
                   (invoke "./configure" (string-append "--prefix=" out))
                   (invoke "make")
                   (invoke "make" "install"))
                 
                 ;; Build PFP
                 (with-directory-excursion "pfp"
                   (invoke "cmake" "-DCMAKE_INSTALL_PREFIX=${out}" ".")
                   (invoke "make")
                   (invoke "make" "install"))
                 
                 ;; Build BigRePair
                 (with-directory-excursion "bigrepair"
                   (invoke "make")
                   (copy-file "procdic" (string-append out "/bin/procdic"))
                   (copy-file "postproc" (string-append out "/bin/postproc"))
                   (copy-file "largeb_repair/largeb_irepair" (string-append out "/bin/largeb_irepair")))
                 
                 ;; Build PFP-Thresholds
                 (with-directory-excursion "pfp-thresholds"
                   (invoke "cmake" "-DCMAKE_INSTALL_PREFIX=${out}" ".")
                   (invoke "make")
                   (invoke "make" "install")))
               #t)))
         (add-after 'install 'wrap-python-script
           (lambda* (#:key inputs outputs #:allow-other-keys)
             (let* ((out (assoc-ref outputs "out"))
                    (python (assoc-ref inputs "python"))
                    (script (string-append out "/bin/moni")))
               (wrap-program script
                 `("PYTHONPATH" ":" prefix (,(string-append out "/lib/python"
                                                            ,(version-major+minor (package-version python))
                                                            "/site-packages")))
                 `("PATH" ":" prefix (,(string-append out "/bin"))))
               #t))))))
    (inputs
     (list zlib
           bzip2
           xz
           gnutls
           openssl
           ncurses
           python-3
           curl
           gcc-12
           autoconf
           automake
           boost
           nss-certs))
    (native-inputs
     (list cmake
           git
           pkg-config))
    (synopsis "Moni-Align")
    (description "Moni-Align project")
    (home-page "https://github.com/maxrossi91/moni-align")
    (license license:expat)))

moni-align-git
