import React from "react";

import Base from "../../components/base"
import {Link} from "gatsby";


export default function Index() {
    return (
        <Base>
            <h1 className="">This is the guides page</h1>
            <p>This is an example <mark>paragraph</mark></p>

            <div className="list-group">
                <Link to="/guides" className="d-flex justify-content-between align-items-center list-group-item list-group-item-action p-3">
                    <span className="">How to play E-major</span>
                    <div>
                        <i className="bi bi-star-fill"></i>
                        <i className="bi bi-star-fill"></i>
                        <i className="bi bi-star-fill"></i>
                        <i className="bi bi-star-fill"></i>
                        <i className="bi bi-star-half"></i>
                    </div>
                </Link>
                <Link to="#" className="list-group-item list-group-item-action">A third link item</Link>
                <Link to="#" className="list-group-item list-group-item-action">A fourth link item</Link>
                <Link to="#" className="list-group-item list-group-item-action disabled">A disabled link item</Link>
            </div>

        </Base>
    );
}