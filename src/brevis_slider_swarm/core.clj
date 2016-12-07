(ns brevis-slider-swarm.core
  (:gen-class)
  (:require [seesaw.core :as seesaw]
            [seesaw.mig :as mig])
  (:use [brevis.graphics.basic-3D]
        [brevis.physics collision core space utils]
        [brevis.shape box sphere cone]
        [brevis core osd vector camera utils display image random]))

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

(def num-birds (atom 500))

(def centering-weight (atom 10))
(def avoidance-weight (atom 10))
(def straying-weight (atom 1))

(def boundary (atom 300))

(def max-velocity (atom 5))
(def max-acceleration (atom 10))

(def control-frame (atom nil))

(defn make-frame []
  (seesaw/frame 
    :title "Brevis Swarm"
    :content
    (mig/mig-panel
      :items [["<html>Slide the sliders to change the weights for the swarm</html>" "span, growx"]
              ["Centering" "gap 10"]
              [(seesaw/slider :id :centering   :min -500 :max 500 :paint-ticks? true :major-tick-spacing 200 :paint-labels? true) "span, growx"]
              ["Avoidance" "gap 10"] 
              [(seesaw/slider :id :avoidance :min -500 :max 500 :paint-ticks? true :major-tick-spacing 200 :paint-labels? true) "span, growx"]
              ["Straying" "gap 10"]
              [(seesaw/slider :id :straying  :min -500 :max 500 :paint-ticks? true :major-tick-spacing 200 :paint-labels? true) "span, growx"]])))

(defn update-weights [root]
  (let [{:keys [centering avoidance straying]} (seesaw/value root)] ; <- Use (value) to get map of values
    (reset! centering-weight centering)
    (reset! avoidance-weight avoidance)
    (reset! straying-weight straying)))

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
  (move (make-real {:type :bird
                    :color (vec4 1 0 0 1)
                    :shape (create-cone 10.2 1.5)})
        position))
  
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

(defn fly
  "Change the acceleration of a bird."
  [bird]
  (let [nbrs (filter bird? (get-neighbor-objects bird))
        bird-pos (get-position bird)

        neighbor-positions (map get-position nbrs)
        neighbor-center (if (zero? (count neighbor-positions))
                          (vec3 0 0 0)
                          (mul-vec3 (apply add-vec3 neighbor-positions)
                                    (/ (count neighbor-positions))))
        
        closest-bird (get-closest-neighbor bird)
        
        new-acceleration (add-vec3 (mul-vec3 (vec3 (- (lrand) 0.5) (- (lrand) 0.5) (- (lrand) 0.5))
                                             @straying-weight)
                                   (mul-vec3 (sub-vec3 bird-pos neighbor-center)
                                             @centering-weight)
                                   (if closest-bird
                                     (mul-vec3 (sub-vec3 bird-pos (get-position closest-bird))
                                               @avoidance-weight)
                                     (vec3 0 0 0)))]    
    (set-velocity
      (set-acceleration
        (if (or (> (java.lang.Math/abs (x-val bird-pos)) @boundary) 
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

(add-collision-handler :bird :bird bump)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ## brevis control code

(defn initialize-simulation
  "This is the function where you add everything to the world."
  []  
  (init-world)
  (init-view)  

  (set-camera-information (vec3 -10.0 57.939613 -890.0) (vec4 1.0 0.0 0.0 0.0))
  
  (let [root (make-frame)]
    (seesaw/listen (map #(seesaw/select root [%]) [:#centering :#avoidance :#straying]) :change
                   (fn [e]
                     (update-weights root)))
    (reset! control-frame root)
    (seesaw/invoke-later
      (-> root
        seesaw/pack!
        seesaw/show!)))
  
  (add-destroy-hook (fn [] (seesaw/dispose! @control-frame))) 
  
  (set-dt 1)
  (set-neighborhood-radius 50)
  (default-display-text)
  (dotimes [_ @num-birds]
    (add-object (random-bird))))

;; Start zee macheen
(defn -main [& args]
  (if-not (empty? args)
    (start-nogui initialize-simulation)
    (start-gui initialize-simulation)))

(autostart-in-repl -main)
