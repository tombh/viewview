#[tokio::main]
async fn main() {
    let static_directory = tower_http::services::ServeDir::new("./static");
    let serve_index = tower_http::services::ServeFile::new("./scripts/tile_packing.html");

    let app = axum::Router::new()
        .fallback_service(serve_index)
        .nest_service("/static", static_directory)
        .layer(tower::ServiceBuilder::new().layer(tower_http::trace::TraceLayer::new_for_http()));

    // run our app with hyper, listening globally on port 3000
    let listener = tokio::net::TcpListener::bind("0.0.0.0:3333").await.unwrap();
    axum::serve(listener, app).await.unwrap();
}
