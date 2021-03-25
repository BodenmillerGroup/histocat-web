import React from "react";

type StatusProps = {
  info?: any;
  warn?: any;
  removeGridComponent: any;
};

export default function Status(props: StatusProps) {
  const warnClass = "alert alert-warning my-0 details";
  const { info, warn } = props;
  const messages: JSX.Element[] = [];
  if (info) {
    messages.push(
      <p className="details" key="info">
        {info}
      </p>
    );
  }
  if (warn) {
    messages.push(
      <p className={warnClass} key="warn">
        {warn}
      </p>
    );
  }
  return <React.Fragment>{messages}</React.Fragment>;
}
