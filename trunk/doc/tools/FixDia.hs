module Main where

import Control.Arrow
import Data.List
import Data.String
import System
import Text.Regex

import Debug.Trace

t x = trace (show x) x `seq` x

subs = lines >>> split [""] >>> map t >>> map sub
  where sub [pat,rep] inp = subRegex (mkRegex pat) inp rep

xform subs = lines >>> map f >>> unlines
  where f line =
          let parts = split "\\$" line
              xformed [a,b,c] = [a, change b, c]
          in if length parts == 3
               then join "$" (xformed parts)
               else line
        change x  = foldl' apply x subs
        apply x f = f x

main = do
  args  <- getArgs
  rules <- readFile $ args !! 0
  interact $ xform $ subs rules
