import React from "react";

type TooltipContentProps = {
  info: object;
};

export default function TooltipContent(props: TooltipContentProps) {
  const { info } = props;

  return (
    <table>
      <tbody>
        {Object.entries(info).map(([key, value]) => (
          <tr key={key}>
            <th>{key}</th>
            <td>{value}</td>
          </tr>
        ))}
      </tbody>
    </table>
  );
}
