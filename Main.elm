port module Main exposing (..)

import Browser
import Html exposing (Html, button, div, input, span, text)
import Html.Attributes exposing (disabled, placeholder, style, type_, value)
import Html.Events exposing (onClick, onInput)
import Json.Decode as D
import Time


-- PORTS


{-| Send a number string to Rust; Rust replies on queryResult. -}
port query : String -> Cmd msg


port queryResult : (String -> msg) -> Sub msg


{-| Ask JS to find the next/previous prime. JS replies on nextPrimeFound / prevPrimeFound. -}
port findNextPrime : String -> Cmd msg


port nextPrimeFound : (String -> msg) -> Sub msg


port findPrevPrime : String -> Cmd msg


port prevPrimeFound : (String -> msg) -> Sub msg



-- TYPES


type alias Factor =
    { value : String
    , exp : Int
    , isPrime : Bool
    }


type QueryResult
    = InputError String
    | QueryOk
        { primality : String -- "prime" | "composite" | "probably_prime" | "checking"
        , probability : Float
        , factorsComplete : Bool
        , factors : List Factor
        }


-- JSON DECODER


factorDecoder : D.Decoder Factor
factorDecoder =
    D.map3 Factor
        (D.field "value" D.string)
        (D.field "exp" D.int)
        (D.field "is_prime" D.bool)


queryResultDecoder : D.Decoder QueryResult
queryResultDecoder =
    D.oneOf
        [ D.map InputError (D.field "error" D.string)
        , D.map4
            (\prim prob complete facs ->
                QueryOk
                    { primality = prim
                    , probability = prob
                    , factorsComplete = complete
                    , factors = facs
                    }
            )
            (D.field "primality" D.string)
            (D.field "probability" D.float)
            (D.field "factors_complete" D.bool)
            (D.field "factors" (D.list factorDecoder))
        ]


parseResult : String -> QueryResult
parseResult json =
    case D.decodeString queryResultDecoder json of
        Ok r ->
            r

        Err _ ->
            InputError ("Parse error: " ++ json)



-- MODEL


type alias Model =
    { input : String
    , lastResult : Maybe QueryResult
    , searchingPrime : Bool
    }


init : () -> ( Model, Cmd Msg )
init _ =
    ( { input = "", lastResult = Nothing, searchingPrime = False }
    , Cmd.none
    )



-- UPDATE


type Msg
    = InputChanged String
    | GotQueryResult String
    | GotNextPrime String
    | GotPrevPrime String
    | RequestNextPrime
    | RequestPrevPrime
    | Tick Time.Posix


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    case msg of
        InputChanged newInput ->
            ( { model
                | input = newInput
                , lastResult = Nothing
                , searchingPrime = False
              }
            , query newInput
            )

        GotQueryResult json ->
            ( { model | lastResult = Just (parseResult json) }
            , Cmd.none
            )

        GotNextPrime n ->
            ( { model
                | input = n
                , lastResult = Nothing
                , searchingPrime = False
              }
            , query n
            )

        GotPrevPrime n ->
            ( { model
                | input = n
                , lastResult = Nothing
                , searchingPrime = False
              }
            , query n
            )

        RequestNextPrime ->
            if model.input == "" then
                ( model, Cmd.none )

            else
                ( { model | searchingPrime = True }
                , findNextPrime model.input
                )

        RequestPrevPrime ->
            if model.input == "" then
                ( model, Cmd.none )

            else
                ( { model | searchingPrime = True }
                , findPrevPrime model.input
                )

        Tick _ ->
            if model.input == "" then
                ( model, Cmd.none )

            else
                ( model, query model.input )



-- SUBSCRIPTIONS


subscriptions : Model -> Sub Msg
subscriptions _ =
    Sub.batch
        [ queryResult GotQueryResult
        , nextPrimeFound GotNextPrime
        , prevPrimeFound GotPrevPrime
        , Time.every 100 Tick
        ]



-- VIEW HELPERS


viewPrimality : String -> Float -> Html Msg
viewPrimality primality probability =
    case primality of
        "prime" ->
            span [ style "color" "#006600", style "font-weight" "bold" ]
                [ text "Prime" ]

        "composite" ->
            span [ style "color" "#880000", style "font-weight" "bold" ]
                [ text "Composite" ]

        "probably_prime" ->
            let
                pct =
                    round (probability * 1.0e10) |> toFloat |> (*) 1.0e-8
            in
            span [ style "color" "#004488" ]
                [ text ("Probably prime (confidence: " ++ String.fromFloat pct ++ "%)") ]

        "checking" ->
            span [ style "color" "#888" ] [ text "Checking…" ]

        _ ->
            text primality


viewFactor : Factor -> String
viewFactor f =
    let
        base =
            if f.isPrime then
                f.value

            else
                f.value ++ " [composite]"

        expStr =
            if f.exp == 1 then
                base

            else
                base ++ "^" ++ String.fromInt f.exp
    in
    expStr


viewFactors : Bool -> List Factor -> Html Msg
viewFactors complete factors =
    let
        parts =
            List.map viewFactor factors

        factorStr =
            String.join " × " parts

        label =
            if complete then
                "Factors: "

            else
                "Factoring: "
    in
    div []
        [ text (label ++ factorStr) ]



-- VIEW


view : Model -> Html Msg
view model =
    let
        isValidInput =
            model.input /= ""

        canNavigate =
            isValidInput && not model.searchingPrime
    in
    div
        [ style "font-family" "monospace"
        , style "max-width" "600px"
        , style "margin" "40px auto"
        , style "padding" "20px"
        ]
        [ -- Input row with navigation buttons
          div
            [ style "display" "flex"
            , style "gap" "8px"
            , style "align-items" "center"
            , style "margin-bottom" "12px"
            ]
            [ button
                [ onClick RequestPrevPrime
                , disabled (not canNavigate)
                , style "padding" "6px 12px"
                ]
                [ text "← Prev Prime" ]
            , input
                [ type_ "text"
                , placeholder "Enter a number"
                , value model.input
                , onInput InputChanged
                , style "flex" "1"
                , style "padding" "6px"
                , style "font-family" "monospace"
                , style "font-size" "1em"
                ]
                []
            , button
                [ onClick RequestNextPrime
                , disabled (not canNavigate)
                , style "padding" "6px 12px"
                ]
                [ text "Next Prime →" ]
            ]
        , -- Status display
          case model.lastResult of
            Nothing ->
                if model.searchingPrime then
                    div [ style "color" "#888" ] [ text "Searching for prime…" ]

                else
                    div [] []

            Just (InputError msg) ->
                div [ style "color" "#888" ] [ text msg ]

            Just (QueryOk r) ->
                div []
                    [ div [ style "margin-bottom" "6px" ]
                        [ if model.searchingPrime then
                            span [ style "color" "#888" ] [ text "Searching for prime…" ]

                          else
                            viewPrimality r.primality r.probability
                        ]
                    , viewFactors r.factorsComplete r.factors
                    ]
        ]



-- MAIN


main : Program () Model Msg
main =
    Browser.element
        { init = init
        , view = view
        , update = update
        , subscriptions = subscriptions
        }
