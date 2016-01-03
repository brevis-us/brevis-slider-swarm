(ns brevis-morphoregulation-swarm.core
  (:gen-class)
  (:use [brevis.graphics basic-3D visual-overlays]
        [brevis.physics collision core space utils]
        [brevis.shape box sphere cone]
        [brevis core osd vector camera utils display image parameters plot random]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ## Swarm
;;
;; Swarm simulations are models of flocking behavior in collections of organisms.   
;;
;; For reference see:
;;
;;   Reynolds, Craig W. "Flocks, herds and schools: A distributed behavioral model." ACM SIGGRAPH Computer Graphics. Vol. 21. No. 4. ACM, 1987.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ## Globals

(def num-birds (atom 2000))

(def avoidance-distance (atom 25))
(def boundary (atom 300))

(def max-velocity (atom 5))
(def max-acceleration (atom 10))

(def swarm-center (atom (vec3 0 0 0)))
(def swarm-center-weight (atom 0.0005))
;(def swarm-center-weight (atom 0))
(def swarm-random-weight (atom 0.001))
(def swarm-velocity-matching-weight (atom 0))

(def swarm-far-weight (atom 0.001))
(def swarm-close-weight (atom 10))

; A potential function, for now, takes 2 positions, assumes the first is the focus
;(def swarm-bird-bird-potential-fn (atom (fn [p1 p2]
;                                          (let [l (length-vec3 (sub-vec3 p1 p2))]
;                                            (

(set-param :vegf-max 1
           :DAPT false
           :mutant-state nil
           :notch-norm 1;(* 1000 50); these 50's are scale, might kinda represent a volume/area 
           :dll4-max 1;(* 1000 50)
           :vegfr-max 1
           :vegfr-min 0;(* 6.89 50)
           :k-actvegfr-dll4 0.75; This is Kate's "delta", active-vegfr that gets converted to dll4
           :notch-queue-length 28
           :vegfr-queue-length 28
           :dt 1
           :bird-radius 25)

; https://en.wikipedia.org/wiki/Lennard-Jones_potential, C H_4 values for N_2 (http://www.sandia.gov/~ajasper/pub/lj.pdf)   
;:lj-epsilon 88.7
;:lj-sigma 3.68)
(let [epsilon 18.7
      sigma 3.68
      A (* 4 epsilon (java.lang.Math/pow sigma 12))
      B (* 4 epsilon (java.lang.Math/pow sigma 6))]
  (set-param :lj-epsilon epsilon
             :lj-sigma sigma
             :lj-A A
             :lj-B B))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ## Birds

(defn bird?
  "Is a thing a bird?"
  [thing]
  (= (get-type thing) :bird))

(defn random-bird-position
  "Returns a random valid bird position."
  [] 
  (vec3 (- (rand @boundary) (/ @boundary 2)) 
        (- (rand @boundary) (/ @boundary 2)) 
        (- (rand @boundary) (/ @boundary 2))))

(defn make-bird
  "Make a new bird. At the specified location."
  [position]  
  (assoc (move (make-real {:type :bird
                           :color (vec4 1 0 0 0.5)
                           :shape (create-sphere (get-param :bird-radius)) 
                           #_(create-cone 10.2 1.5)})
               position)
         :notch-queue (vec (repeat (get-param :notch-queue-length) 0))
         :vegfr-queue (vec (repeat (get-param :vegfr-queue-length) 0))))
  
(defn random-bird
  "Make a new random bird."
  []
  (make-bird (random-bird-position)))    

(defn bound-acceleration
  "Keeps the acceleration within a reasonable range."
  [v]  
  (if (> (length v) @max-acceleration)
    (mul (div v (length v)) @max-acceleration)
    v))

(defn bound-velocity
  "Keeps the acceleration within a reasonable range."
  [v]  
  (if (> (length v) @max-velocity)
    (mul (div v (length v)) @max-velocity)
    v))

(defn periodic-boundary
  "Change a position according to periodic boundary conditions."
  [pos]
  (let [x (x-val pos)
        y (y-val pos)
        z (z-val pos)]
    (vec3 (cond (> x @boundary) (- (mod x @boundary) @boundary)
                (< x (- @boundary)) (mod (- x) @boundary)
                :else x)
          (cond (> y @boundary) (- (mod y @boundary) @boundary)
                (< y (- @boundary)) (mod (- y) @boundary)
                :else y)
          (cond (> z @boundary) (- (mod z @boundary) @boundary)
                (< z (- @boundary)) (mod (- z) @boundary)
                :else z))))

(defn fixed-boundary
  "Change a position according to periodic boundary conditions."
  [pos]
  (let [size-adjusted-boundary (- @boundary (get-param :bird-radius))
        x (x-val pos)
        y (y-val pos)
        z (z-val pos)]
    (vec3 (cond (> x size-adjusted-boundary) size-adjusted-boundary
                (< x (- size-adjusted-boundary)) (- size-adjusted-boundary)
                :else x)
          (cond (> y size-adjusted-boundary) size-adjusted-boundary
                (< y (- size-adjusted-boundary)) (- size-adjusted-boundary)
                :else y)
          (cond (> z size-adjusted-boundary) size-adjusted-boundary
                (< z (- size-adjusted-boundary)) (- size-adjusted-boundary)
                :else z))))

(defn fly
  "Change the acceleration of a bird."
  [bird]
  (let [bird-pos (get-position bird)
        
        closest-bird (get-closest-neighbor bird)
        
        ; Radial
        #_env-vegf-current #_(* (- 1 (/ (length-vec3 (sub @swarm-center bird-pos))
                                      @boundary))
                              (get-param :vegf-max))
        
        ; Uniform
        ;env-vegf-current 1
        
        ; Vertical
        env-vegf-current (/ (+ (/ (y-val-vec3 bird-pos) @boundary) 1) 2)
        
        vegfr-tot (min (max (+ #_(first (:notch-queue bird))
                               (let [n 5] (/ (apply + (take n (:notch-queue bird))) n))
                               (lrand 0.0001))
                            (get-param :vegfr-min))
                       (get-param :vegfr-max))
        dll4-tot (if closest-bird
                   (min (+ (* (get-param :k-actvegfr-dll4)
                              (let [n 5] (/ (apply + (take n (:vegfr-queue closest-bird))) n)))
                              ;(first (:vegfr-queue closest-bird)))
                           (lrand 0.0001))
                        (get-param :dll4-max))
                   0)
        
        act-not-current (cond (get-param :DAPT) 0
                              (= (get-param :mutant-state) :perm-inhib) (get-param :notch-norm)
                              (= (get-param :mutant-state) :perm-active) 0
                              :else (- 1 dll4-tot))
        act-vegfr-current (min vegfr-tot env-vegf-current); min means: number bound, whether limited by vegfr or vegf
        
        notch-queue (conj (vec (rest (:notch-queue bird)))
                          act-not-current)
                          ;dll4-tot)
        vegfr-queue (conj (vec (rest (:vegfr-queue bird)))
                          act-vegfr-current)
        
        migration-rate (if closest-bird                         
                         (/ dll4-tot (get-param :dll4-max))
                         1)
        
        new-acceleration (mul (add ;Swarm centering
                                   (mul (sub @swarm-center bird-pos)
                                        @swarm-center-weight)
                                   ;Swarm randomizing
                                   (mul (vec3 (- (rand) 0.5) (- (rand) 0.5) (- (rand) 0.5))
                                        @swarm-random-weight)
                                   ;Velocity matching
                                   (if closest-bird
                                     (mul (sub (get-velocity closest-bird) (get-velocity bird))
                                          @swarm-velocity-matching-weight)
                                     (vec3 0 0 0))
                                   ;Attraction/repulsion
                                   (if closest-bird
                                     ; Lennard-jones: https://en.wikipedia.org/wiki/Lennard-Jones_potential
                                     (let [dvec (sub bird-pos (get-position closest-bird)) 
                                           len (length dvec)
                                           potential (- (/ (get-param :lj-A) (java.lang.Math/pow len 12))
                                                        (/ (get-param :lj-B) (java.lang.Math/pow len 6)))]
                                       (mul-vec3 dvec
                                                 (/ potential len)))
                                     #_(let [dvec (sub bird-pos (get-position closest-bird)) 
                                            len (length dvec)]
                                        (if (<= len @avoidance-distance)
                                            ;; If far from neighbor, get closer
                                            (mul dvec
                                                 @swarm-far-weight)
                                            ;; If too close to neighbor, move away
                                            (mul dvec (- @swarm-close-weight))))
                                     (vec3 0 0 0)))
                              migration-rate);
        new-acceleration (if (zero? (length new-acceleration))
                           new-acceleration
                           (mul new-acceleration (/ migration-rate
                                                    (length new-acceleration))))]
    ;(println new-acceleration)
    (set-velocity
      (set-acceleration
        (set-color (assoc (move bird (fixed-boundary bird-pos))
                          ;(move bird (periodic-boundary bird-pos))
                          ;bird
                          :notch-queue notch-queue
                          :vegfr-queue vegfr-queue)
                   (vec4 (/ act-not-current (get-param :dll4-max)) act-vegfr-current 1 0.5))
        #_(if (or (> (java.lang.Math/abs (x-val bird-pos)) @boundary) 
                 (> (java.lang.Math/abs (y-val bird-pos)) @boundary) 
                 (> (java.lang.Math/abs (z-val bird-pos)) @boundary)) 
           (move bird (periodic-boundary bird-pos) #_(vec3 0 25 0))
           bird)
        (bound-acceleration new-acceleration))
      (bound-velocity (get-velocity bird)))))

(enable-kinematics-update :bird); This tells the simulator to move our objects
(add-update-handler :bird fly); This tells the simulator how to update these objects

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ## Collision handling
;;
;; Collision functions take [collider collidee] and return [collider collidee]
;; This is only called once per pair of colliding objects.

(defn bump
  "Collision between two birds."
  [bird1 bird2]  
  [(set-color bird1 (vec4 (rand) (rand) (rand) 1))
   bird2])

#_(add-collision-handler :bird :bird bump)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ## brevis control code

(defn draw-bounding-box
  "Draw a bounding box as a brevis visual overlay."
  [xmin xmax ymin ymax zmin zmax]
  (let [corners (for [x [xmin xmax] y [ymin ymax] z [zmin zmax]] (vec3 x y z))]    
    (loop [sources corners]
      (when (> (count sources) 1)
        (doseq [dest (rest sources)]
          (add-line (first sources) dest))
        (recur (rest sources))))))

(defn initialize-simulation
  "This is the function where you add everything to the world."
  []  
  (init-world)
  (init-view)  

  (set-camera-information (vec3 -10.0 57.939613 -890.0) (vec4 1.0 0.0 0.0 0.0))
  
  (add-plot-handler
   (fn [] 
     (let [birds (filter bird? (all-objects))]
       [(* (get-time) (get-dt)) 
        (/ (apply + (map #(first (:notch-queue %)) birds))
           (count birds))]))
   :interval 50
   :title "Average notch")
  
  (add-plot-handler
   (fn [] 
     (let [birds (filter bird? (all-objects))]
       [(* (get-time) (get-dt)) 
        (/ (apply + (map #(first (:vegfr-queue %)) birds))
           (count birds))]))
   :interval 50
   :title "Average vegfr")
  
  (set-dt (get-param :dt))
  (set-neighborhood-radius 50)
  (default-display-text)
  (set-parallel true)
  (draw-bounding-box (- @boundary) @boundary  (- @boundary) @boundary  (- @boundary) @boundary)
  
  (dotimes [_ @num-birds]
    (add-object (random-bird))))

;; Start zee macheen
(defn -main [& args]
  (if-not (empty? args)
    (start-nogui initialize-simulation)
    (start-gui initialize-simulation)))

(autostart-in-repl -main)
