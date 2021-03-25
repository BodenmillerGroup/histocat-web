import React from "react";
import expect from "expect";
import VitessceGridLayout from "./VitessceGridLayout";
import { render, screen } from "@testing-library/react";

describe("VitessceGridLayout.tsx", () => {
  describe("<VitessceGridLayout />", () => {
    function FakeComponent(props: any) {
      const { text } = props;
      return <span>{text}</span>;
    }
    /* eslint-disable object-curly-newline */
    /* eslint-disable object-property-newline */
    const layoutJson = {
      columns: {
        600: [0, 2, 4, 8],
      },
      components: [{ component: "FakeComponent", props: { text: "Hello World" }, x: 0, y: 0, w: 2 }],
    };
    /* eslint-enable */
    it("mount() works", () => {
      const { container } = render(
        <VitessceGridLayout layout={layoutJson} getComponent={() => FakeComponent} draggableHandle=".my-handle" />
      );

      expect(screen.getAllByText("Hello World")).toHaveLength(1);
      // expect(wrapper.find('.react-grid-item').length).toEqual(1);
      // expect(wrapper.find('.react-grid-item').text()).toEqual('Hello World');
      // expect(wrapper.find('.react-grid-item span:not(.react-resizable-handle)').length).toEqual(1);

      // const style = wrapper.find('style');
      // expect(style.length).toEqual(1);
      // expect(style.text()).toContain('.my-handle {');
      // expect(style.text()).toContain('.my-handle:active {');
    });

    it("rowHeight works", () => {
      const { container } = render(
        <VitessceGridLayout
          layout={layoutJson}
          getComponent={() => FakeComponent}
          draggableHandle=".my-handle"
          rowHeight={123}
        />
      );

      expect(container.getElementsByClassName("react-grid-item")[0].style.height).toEqual("123px");
    });
  });
});
